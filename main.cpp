/*****************************************************************************
 * VCLib                                                                     *
 * Visual Computing Library                                                  *
 *                                                                           *
 * Copyright(C) 2021-2026                                                    *
 * Visual Computing Lab                                                      *
 * ISTI - Italian National Research Council                                  *
 *                                                                           *
 * All rights reserved.                                                      *
 *                                                                           *
 * This program is free software; you can redistribute it and/or modify      *
 * it under the terms of the Mozilla Public License Version 2.0 as published *
 * by the Mozilla Foundation; either version 2 of the License, or            *
 * (at your option) any later version.                                       *
 *                                                                           *
 * This program is distributed in the hope that it will be useful,           *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              *
 * Mozilla Public License Version 2.0                                        *
 * (https://www.mozilla.org/en-US/MPL/2.0/) for more details.                *
 ****************************************************************************/

#include <vclib/embree/scene.h>
#include <vclib/igl/booleans.h>
#include <vclib/io.h>
#include <vclib/meshes.h>
#include <vclib/algorithms/core/fibonacci.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <limits>
#include <map>
#include <numeric>
#include <utility>
#include <vector>

#include "include/struct.h"
#include "include/helper.h"
#include "include/functions.h"
#include "include/debug_output.h"

struct MoldCheckMetrics
{
    double score = -std::numeric_limits<double>::infinity();
    double hitRatio = 0.0;
    double compactness = 0.0;
    vcl::uint hitCount = 0;
    double reduceRatio = 0.0;
    double hiddenRatio = 0.0;
};

static double moldQualityScore(
    double hitRatio,
    double compactness,
    double reduceRatio,
    double hiddenRatio)
{
    return
        0.35 * hitRatio +
        0.25 * compactness +
        0.20 * (1 - (hiddenRatio)) +
        0.20 * (1 - (reduceRatio));
}

static vcl::PolyMesh makeNormalSideCutterVolume(
    const vcl::TriMesh& moldSurfaceMesh,
    vcl::Point3d        normalDirection,
    double              extrusionDistance)
{
    using namespace vcl;

    normalDirection.normalize();

    PolyMesh cutterVolume;

    std::vector<uint> bottomIds;
    std::vector<uint> topIds;

    bottomIds.reserve(moldSurfaceMesh.vertexCount());
    topIds.reserve(moldSurfaceMesh.vertexCount());

    for (const auto& v : moldSurfaceMesh.vertices()) {
        const Point3d p = v.position();

        bottomIds.push_back(cutterVolume.addVertex(p));
        topIds.push_back(cutterVolume.addVertex(
            p + normalDirection * extrusionDistance));
    }

    using EdgeKey = std::pair<uint, uint>;

    struct BoundaryEdgeData
    {
        uint a = 0;
        uint b = 0;
        uint count = 0;
    };

    std::map<EdgeKey, BoundaryEdgeData> edgeMap;

    auto addBoundaryCandidateEdge = [&](uint a, uint b) {
        const EdgeKey key = std::minmax(a, b);
        auto& data = edgeMap[key];

        if (data.count == 0) {
            data.a = a;
            data.b = b;
        }

        ++data.count;
    };

    for (const auto& f : moldSurfaceMesh.faces()) {
        uint a = f.vertex(0)->index();
        uint b = f.vertex(1)->index();
        uint c = f.vertex(2)->index();

        const Point3d pa = moldSurfaceMesh.vertex(a).position();
        const Point3d pb = moldSurfaceMesh.vertex(b).position();
        const Point3d pc = moldSurfaceMesh.vertex(c).position();

        Point3d faceNormal = (pb - pa).cross(pc - pa);

        if (faceNormal.dot(normalDirection) < 0.0) {
            std::swap(b, c);
        }

        const uint ba = bottomIds[a];
        const uint bb = bottomIds[b];
        const uint bc = bottomIds[c];

        const uint ta = topIds[a];
        const uint tb = topIds[b];
        const uint tc = topIds[c];

        cutterVolume.addFace(ta, tb, tc);

        cutterVolume.addFace(bc, bb, ba);

        addBoundaryCandidateEdge(a, b);
        addBoundaryCandidateEdge(b, c);
        addBoundaryCandidateEdge(c, a);
    }

    for (const auto& [_, edge] : edgeMap) {
        if (edge.count != 1) {
            continue;
        }

        const uint a = edge.a;
        const uint b = edge.b;

        const uint ba = bottomIds[a];
        const uint bb = bottomIds[b];

        const uint ta = topIds[a];
        const uint tb = topIds[b];

        cutterVolume.addFace(ba, bb, tb);
        cutterVolume.addFace(ba, tb, ta);
    }

    updateBoundingBox(cutterVolume);

    return cutterVolume;
}

MoldCheckMetrics moldCheck(
    vcl::PolyMesh              m,
    const std::vector<double>& gridCellSideLengths,
    bool                       debug,
    vcl::Point3d               direction,
    const double               coneAngleDegrees,
    const double               draftAngleDegrees,
    const double               marginFactor,
    const std::string&         debugResultsSubdir = "",
    const vcl::PolyMesh*       moldMesh = nullptr)
{
    using namespace vcl;

    const double CONE_COS_THRESHOLD =
        std::cos(coneAngleDegrees * M_PI / 180.0);

    if (debug) {
        std::cout << "=== moldCheck started ===\n";
        std::cout.flush();
    }

    updateBoundingBox(m);

    const double MAX_DISTANCE = m.boundingBox().diagonal();
    const float EPS = 1e-12f * MAX_DISTANCE;
    const float RAY_EPS = 1e-6f * MAX_DISTANCE;

    embree::Scene scene(m);

    direction.normalize();

    double minProj = std::numeric_limits<double>::infinity();
    for (const auto& vv : m.vertices()) {
        minProj = std::min(minProj, vv.position().dot(direction));
    }

    const Point3d planePoint = direction * minProj;
    const Planed plane(planePoint, direction);
    const double margin = marginFactor * MAX_DISTANCE;

    GridChoice grid;

    const auto [u, v] = makePlane(
        m,
        plane,
        planePoint,
        direction,
        margin,
        EPS,
        grid);

    const double lenU = grid.maxU - grid.minU;
    const double lenV = grid.maxV - grid.minV;

    if (lenU <= EPS || lenV <= EPS) {
        return {};
    }

    makeGrid(grid, gridCellSideLengths);

    const double cellDu = grid.sideU;
    const double cellDv = grid.sideV;
    const double cellArea = cellDu * cellDv;

    std::vector<uint> allCells(grid.rows * grid.cols);
    std::iota(allCells.begin(), allCells.end(), 0);

    std::vector<CellData> cells(allCells.size());

    parallelFor(allCells, [&](uint idx) {
        const CellData cell = makeCellGeometry(idx, grid, planePoint, u, v);
        cells[idx] = shootRayOnCell(
            cell,
            m,
            scene,
            planePoint,
            direction,
            MAX_DISTANCE,
            RAY_EPS);
    });

    uint rawHitCount = 0;
    std::vector<uint> hitCellIds;
    hitCellIds.reserve(cells.size());

    for (uint i = 0; i < cells.size(); ++i) {
        if (cells[i].hasHit) {
            ++rawHitCount;
            hitCellIds.push_back(i);
        }
    }

    if (debug) {

        std::cout << "Ray casting complete. Hit cells: "
                  << rawHitCount << "/" << allCells.size() << "\n";
        std::cout << "Beginning Clamping phase...\n";
        std::cout.flush();
    }

    parallelFor(allCells, [&](uint idx) {
        computeClampedCell(
            idx,
            cells,
            hitCellIds,
            planePoint,
            direction,
            CONE_COS_THRESHOLD,
            EPS);
    });

    const double REDUCE_POINTS_REFERENCE_CELL_SIDE = 0.4;
    const double reducePointsCellScale =
        ((grid.sideU + grid.sideV) * 0.5) /
        REDUCE_POINTS_REFERENCE_CELL_SIDE;

    const double REDUCE_POINTS_DISTANCE_THRESHOLD =
        0.02 * MAX_DISTANCE * reducePointsCellScale;

    cells = reducePoints(
        cells,
        grid,
        direction,
        draftAngleDegrees,
        EPS,
        REDUCE_POINTS_DISTANCE_THRESHOLD,
        MAX_DISTANCE);

    double totalAreaHit = 0.0;
    double clampedAreaHit = 0.0;
    double hiddenAreaHit = 0.0;
    uint reducedHitCount = 0;

    for (uint i = 0; i < cells.size(); ++i) {
        if (cells[i].hasHit) {
            ++reducedHitCount;
            totalAreaHit += cellArea;
        }

        if (cells[i].hasHit &&
            cells[i].hasClampedHit &&
            std::abs(cells[i].clampedDistance - cells[i].distance) > EPS) {
            clampedAreaHit += cellArea;
        }

        if (cells[i].hasHit && cells[i].hitPoints.size() > 2) {
            hiddenAreaHit += cellArea;
        }
    }

    const double percentClamped =
        (totalAreaHit > 0.0) ?
            (clampedAreaHit / totalAreaHit):
            0.0;

    const double reduceRatio = 
        (rawHitCount > 0) ?
            (reducedHitCount / static_cast<double>(rawHitCount)):
            0.0;

    const double hiddenRatio =
        (totalAreaHit > 0.0) ?
            (hiddenAreaHit / totalAreaHit):
            0.0;

    const HitCellShapeData hitShape = hitCellShape(cells, grid);

    const double hitRatio =
        (cells.size() > 0) ?
            static_cast<double>(reducedHitCount) / cells.size() :
            0.0;

    const MoldCheckMetrics metrics{
        moldQualityScore(
            hitRatio,
            hitShape.compactness,
            reduceRatio,
            hiddenRatio),
        hitRatio,
        hitShape.compactness,
        reducedHitCount,
        reduceRatio,
        hiddenRatio};

    if (debug) {
        std::cout << "Clamping and reduction complete. Hit cells after reduction: "
                  << reducedHitCount << "/" << allCells.size() << "\n";
        std::cout.flush();
    }

    std::vector<CellData> depthCells = cells;

    if (debug) {
        std::cout << "Beginning Depth Smoothing phase...\n";

        depthCells =
            makeDepthCells(
                cells,
                direction,
                grid,
                CONE_COS_THRESHOLD,
                EPS,
                debugResultsSubdir);

        std::cout << "Depth smoothing complete.\n";
        std::cout << "Validating clamped cells...\n";
        std::cout.flush();

        PolyMesh violatingPointsMesh =
            validateClampedCells(
                depthCells,
                allCells,
                direction,
                CONE_COS_THRESHOLD,
                EPS);

        PolyMesh hitPointsMesh;
        hitPointsMesh.enablePerVertexColor();

        for (uint i = 0; i < cells.size(); ++i) {
            if (cells[i].distance == MAX_DISTANCE) {
                continue;
            }

            addColoredPoint(
                hitPointsMesh,
                cells[i].cellCenter + direction * cells[i].distance,
                Color::Yellow);
        }

        PolyMesh hitPointsafterReductionMesh;
        hitPointsafterReductionMesh.enablePerVertexColor();

        for (uint i = 0; i < cells.size(); ++i) {
            if (!cells[i].hasHit) {
                continue;
            }

            addColoredPoint(
                hitPointsafterReductionMesh,
                cells[i].cellCenter + direction * cells[i].distance,
                Color::Blue);
        }

        PolyMesh clampedPointsMesh;
        clampedPointsMesh.enablePerVertexColor();

        for (uint i = 0; i < cells.size(); ++i) {
            if (!cells[i].hasHit) {
                continue;
            }

            addColoredPoint(
                clampedPointsMesh,
                cells[i].cellCenter + direction * cells[i].clampedDistance,
                Color::Blue);
        }

        PolyMesh depthPointsMesh;
        depthPointsMesh.enablePerVertexColor();

        for (uint i = 0; i < depthCells.size(); ++i) {
            const Point3d depthPoint =
                depthCells[i].cellCenter +
                direction * depthCells[i].distance;

            Color depthColor = Color::White;

            if (depthCells[i].hasHit) {
                if (depthCells[i].hasClampedHit) {
                    depthColor = Color::Green;
                }
                else {
                    depthColor = Color::Red;
                }
            }

            addColoredPoint(
                depthPointsMesh,
                depthPoint,
                depthColor);
        }

        const TriMesh planeMesh =
            makeDebugPlaneMesh(grid, planePoint, u, v);

        const TriMesh moldSurfaceMesh =
            createMoldSurface(depthCells, grid, direction);

        PolyMesh moldSurfaceCutterVolume;
        PolyMesh moldNormalSidePiece;

        bool hasMoldNormalSidePiece = false;

        if (moldMesh != nullptr && moldSurfaceMesh.faceCount() > 0) {
            std::cout << "Cutting moldMesh with mold_surface normal side...\n";
            std::cout.flush();

            const double cutterExtrusionDistance = 2.0 * MAX_DISTANCE;

            moldSurfaceCutterVolume =
                makeNormalSideCutterVolume(
                    moldSurfaceMesh,
                    -direction,
                    cutterExtrusionDistance);

            if (moldSurfaceCutterVolume.faceCount() > 0) {
                moldNormalSidePiece =
                    vcl::igl::meshBoolean(
                        *moldMesh,
                        moldSurfaceCutterVolume,
                        vcl::igl::MeshBoolean::INTERSECTION);

                hasMoldNormalSidePiece = true;
            }
        }

        const std::filesystem::path debugOutputDir =
            std::filesystem::path(RESULTS_PATH) /
            debugResultsSubdir;

        std::filesystem::create_directories(debugOutputDir);

        const std::string base =
            (debugOutputDir / "mold_check").string();

        saveMesh(hitPointsMesh, base + "_hit_points.ply");
        saveMesh(
            hitPointsafterReductionMesh,
            base + "_hit_points_after_reduction.ply");
        saveMesh(clampedPointsMesh, base + "_lipschitz_points.ply");
        saveMesh(depthPointsMesh, base + "_mold_points.ply");
        saveMesh(planeMesh, base + "_plane_tested.ply");
        saveMesh(moldSurfaceMesh, base + "_mold_surface.ply");
        saveMesh(violatingPointsMesh, base + "_non-lipschitz_points.ply");

        if (moldSurfaceCutterVolume.faceCount() > 0) {
            saveMesh(
                moldSurfaceCutterVolume,
                base + "_mold_surface_cutter_volume.ply");
        }

        if (hasMoldNormalSidePiece) {
            saveMesh(
                moldNormalSidePiece,
                base + "_mold_piece_1.ply");
        }

        std::cout << "Clamped points: "
                  << clampedPointsMesh.vertexCount() << "\n";
        std::cout << "Depth points: "
                  << depthPointsMesh.vertexCount() << "\n";
        std::cout << "Mold surface points: "
                  << moldSurfaceMesh.vertexCount() << "\n";
        std::cout << "Mold surface faces: "
                  << moldSurfaceMesh.faceCount() << "\n";

        if (moldMesh != nullptr) {
            std::cout << "Mold normal-side piece vertices: "
                      << moldNormalSidePiece.vertexCount() << "\n";
            std::cout << "Mold normal-side piece faces: "
                      << moldNormalSidePiece.faceCount() << "\n";
        }

        std::cout << "Hit cells area: "
                  << hitShape.area << "\n";
        std::cout << "Hit cells perimeter: "
                  << hitShape.perimeter << "\n";
        std::cout << "Hit cells compactness: "
                  << hitShape.compactness << "\n";
        std::cout << "TotalAreaHit: "
                  << totalAreaHit << "\n";
        std::cout << "RawHitCount: "
                  << rawHitCount << "\n";
        std::cout << "ClampedAreaHit: "
                  << clampedAreaHit << "\n";
        std::cout << "percentClamped: "
                  << percentClamped << "\n";
        std::cout << "hiddenAreaHit: "
                  << hiddenAreaHit << "\n";
        std::cout << "hiddenRatio: "
                  << hiddenRatio << "\n";
        std::cout << "reduceRatio: "
                  << reduceRatio << "\n";
        std::cout << "hitRatio: "
                  << hitRatio << "\n";
        std::cout << "hitCount: "
                  << reducedHitCount << "\n";
        std::cout << "qualityScore: "
                  << metrics.score << "\n";

        std::cout << "Saved debug meshes:\n"
                  << " - " << base << "_hit_points.ply\n"
                  << " - " << base << "_hit_points_after_reduction.ply\n"
                  << " - " << base << "_lipschitz_points.ply\n"
                  << " - " << base << "_mold_points.ply\n"
                  << " - " << base << "_plane_tested.ply\n"
                  << " - " << base << "_mold_surface.ply\n"
                  << " - " << base << "_mold_surface_cutter_volume.ply\n"
                  << " - " << base << "_non-lipschitz_points.ply\n";

        if (hasMoldNormalSidePiece) {
            std::cout << " - " << base
                      << "_mold_piece_1.ply\n";
        }

        std::cout << "=== moldCheck completed successfully ===\n";
        std::cout.flush();
    }

    return metrics;
}

static vcl::PolyMesh makeContainingBoxMesh(
    vcl::PolyMesh mesh,
    double marginFactor)
{
    using namespace vcl;
    using namespace vcl::igl;

    updateBoundingBox(mesh);

    const double maxDistance = mesh.boundingBox().diagonal();
    const double boxMargin = marginFactor * maxDistance;

    const Point3d boxMin =
        mesh.boundingBox().min() -
        Point3d(boxMargin, boxMargin, boxMargin);

    const Point3d boxMax =
        mesh.boundingBox().max() +
        Point3d(boxMargin, boxMargin, boxMargin);

    PolyMesh boxMesh;

    const uint v0 =
        boxMesh.addVertex(Point3d(boxMin.x(), boxMin.y(), boxMin.z()));
    const uint v1 =
        boxMesh.addVertex(Point3d(boxMax.x(), boxMin.y(), boxMin.z()));
    const uint v2 =
        boxMesh.addVertex(Point3d(boxMax.x(), boxMax.y(), boxMin.z()));
    const uint v3 =
        boxMesh.addVertex(Point3d(boxMin.x(), boxMax.y(), boxMin.z()));

    const uint v4 =
        boxMesh.addVertex(Point3d(boxMin.x(), boxMin.y(), boxMax.z()));
    const uint v5 =
        boxMesh.addVertex(Point3d(boxMax.x(), boxMin.y(), boxMax.z()));
    const uint v6 =
        boxMesh.addVertex(Point3d(boxMax.x(), boxMax.y(), boxMax.z()));
    const uint v7 =
        boxMesh.addVertex(Point3d(boxMin.x(), boxMax.y(), boxMax.z()));

    boxMesh.addFace(v0, v3, v2);
    boxMesh.addFace(v0, v2, v1);

    boxMesh.addFace(v4, v5, v6);
    boxMesh.addFace(v4, v6, v7);

    boxMesh.addFace(v0, v1, v5);
    boxMesh.addFace(v0, v5, v4);

    boxMesh.addFace(v1, v2, v6);
    boxMesh.addFace(v1, v6, v5);

    boxMesh.addFace(v2, v3, v7);
    boxMesh.addFace(v2, v7, v6);

    boxMesh.addFace(v3, v0, v4);
    boxMesh.addFace(v3, v4, v7);

    updateBoundingBox(boxMesh);

    PolyMesh result =
        meshBoolean(
            boxMesh,
            mesh,
            MeshBoolean::DIFFERENCE);

    return result;
}

int main()
{
    using namespace vcl;

    const auto startTime = std::chrono::steady_clock::now();

    std::filesystem::create_directories(RESULTS_PATH);

    //bigger than boxMesh for it to completely encapsulate the mold mesh
	const double marginFactor = 0.2;
    const uint NUM_PLANES = 100;

    std::vector<Point3d> fibNormals =
        sphericalFibonacciPointSet<Point3d>(NUM_PLANES);

    PolyMesh m =
        loadMesh<PolyMesh>(MESHES_PATH "/bimba_enlarged.ply");

    const PolyMesh moldMesh =
        makeContainingBoxMesh(m, 0.05);

    saveMesh(moldMesh, RESULTS_PATH "/mold-bool.ply");

    std::vector<double> gridCellSideLengths = {0.4, 0.4};

    const double coneAngleDegrees = 5.0;
    const double draftAngleDegrees = 20.0;

    MoldCheckMetrics result;
    MoldCheckMetrics bestResult;
    MoldCheckMetrics worstResult;

    worstResult.score = std::numeric_limits<double>::infinity();

    int bestDirectionIndex = 0;
    int worstDirectionIndex = 0;

    std::vector<std::pair<double, int>> scoredDirections;

    for (const auto& direction : fibNormals) {
        const int directionIndex =
            static_cast<int>(&direction - &fibNormals[0]);

        std::cout << "Processing direction: "
                  << direction << "\n";

        result =
            moldCheck(
                m,
                gridCellSideLengths,
                false,
                direction,
                coneAngleDegrees,
                draftAngleDegrees,
                marginFactor);

        scoredDirections.push_back(
            {result.score, directionIndex});

        std::cout << "Score: " << result.score
                  << " (hit ratio: " << result.hitRatio
                  << ", compactness: " << result.compactness
                  << ", hits: " << result.hitCount
                  << ", reduced: " << result.reduceRatio
                  << ", hidden: " << result.hiddenRatio << ")\n";

        if (result.score > bestResult.score) {
            bestResult = result;
            bestDirectionIndex = directionIndex;

            std::cout << "New best direction found! Index: "
                      << bestDirectionIndex
                      << ", Score: "
                      << bestResult.score
                      << "\n";
        }

        if (result.score < worstResult.score) {
            worstResult = result;
            worstDirectionIndex = directionIndex;

            std::cout << "New worst direction found! Index: "
                      << worstDirectionIndex
                      << ", Score: "
                      << worstResult.score
                      << "\n";
        }
    }

    std::sort(
        scoredDirections.begin(),
        scoredDirections.end(),
        [](const auto& a, const auto& b) {
            return a.first < b.first;
        });

    if (!fibNormals.empty()) {
        std::cout << "Processing debug direction 0\n";
        result =
            moldCheck(
                m,
                gridCellSideLengths,
                true,
                fibNormals[0],
                coneAngleDegrees,
                draftAngleDegrees,
                marginFactor,
                "direction_0",
                &moldMesh);
    }

    result =
        moldCheck(
            m,
            gridCellSideLengths,
            true,
            fibNormals[bestDirectionIndex],
            coneAngleDegrees,
            draftAngleDegrees,
            marginFactor,
            "best",
            &moldMesh);

    result =
        moldCheck(
            m,
            gridCellSideLengths,
            true,
            fibNormals[worstDirectionIndex],
            coneAngleDegrees,
            draftAngleDegrees,
            marginFactor,
            "worst",
            &moldMesh);

    const int medianDebugCount =
        std::min<int>(
            15,
            static_cast<int>(scoredDirections.size()));

    const int medianCenter =
        static_cast<int>(scoredDirections.size() / 2);

    const int medianStart =
        std::max(
            0,
            std::min(
                medianCenter - medianDebugCount / 2,
                static_cast<int>(scoredDirections.size()) -
                    medianDebugCount));

    for (int i = 0; i < medianDebugCount; ++i) {
        const auto& [medianScore, medianDirectionIndex] =
            scoredDirections[medianStart + i];

        const std::string medianDir =
            "median " + std::to_string(i + 1);

        std::cout << "Processing median direction "
                  << (i + 1)
                  << "/"
                  << medianDebugCount
                  << ". Index: "
                  << medianDirectionIndex
                  << ", Score: "
                  << medianScore
                  << ", Output: "
                  << medianDir
                  << "\n";

        result =
            moldCheck(
                m,
                gridCellSideLengths,
                true,
                fibNormals[medianDirectionIndex],
                coneAngleDegrees,
                draftAngleDegrees,
                marginFactor,
                medianDir,
                &moldMesh);
    }

    const auto endTime =
        std::chrono::steady_clock::now();

    const auto elapsedMs =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            endTime - startTime);

    std::cout << "moldCheck execution time: "
              << elapsedMs.count()
              << " ms\n";

    std::cout.flush();

    return 0;
}
