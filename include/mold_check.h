#ifndef MOLD_CHECK_MOLD_CHECK_H
#define MOLD_CHECK_MOLD_CHECK_H

#include <vclib/embree/scene.h>
#include <vclib/algorithms/mesh/update/bounding_box.h>
#include <vclib/igl/booleans.h>
#include <vclib/io.h>
#include <vclib/meshes.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "debug_output.h"
#include "depth.h"
#include "helper.h"
#include "rays.h"
#include "reduce.h"
#include "struct.h"

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

inline MoldCheckMetrics moldCheck(
    vcl::PolyMesh              m,
    const std::vector<double>& gridCellSideLengths,
    bool                       debug,
    vcl::Point3d               direction,
    const double               coneAngleDegrees,
    const double               magentaAngleDegrees,
    const vcl::uint            magentaCellInterval,
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
    const size_t MAX_BIHARMONIC_ITERATIONS = 30000;
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
        cells[idx] = makeCellGeometry(idx, grid, planePoint, u, v);
        cells[idx].distance = MAX_DISTANCE;
		// cells[idx].clampedDistance = MAX_DISTANCE;
        cells[idx].boundaries = {0.0, MAX_DISTANCE};
    });

    if (moldMesh != nullptr) {
        double moldMinProj = std::numeric_limits<double>::infinity();
        for (const auto& vertex : moldMesh->vertices()) {
            moldMinProj = std::min(
                moldMinProj,
                vertex.position().dot(direction));
        }

        const Point3d boundaryPlanePoint = direction * moldMinProj;
        const double boundaryPlaneDistance = minProj - moldMinProj;
        const embree::Scene moldScene(*moldMesh);
        parallelFor(allCells, [&](uint idx) {
            const CellData boundaryCell =
                makeCellGeometry(
                    idx,
                    grid,
                    boundaryPlanePoint,
                    u,
                    v);
            cells[idx].boundaries = shootBoundaryRayOnCell(
                boundaryCell,
                moldScene,
                direction,
                boundaryPlaneDistance,
                MAX_DISTANCE,
                RAY_EPS);
            cells[idx].isDiscarded =
                cells[idx].boundaries[0] == -MAX_DISTANCE &&
                cells[idx].boundaries[1] == MAX_DISTANCE;
        });

        restoreBoundaryCollarCells(
            cells,
            allCells,
            grid,
            MAX_DISTANCE,
            1);

    }

    const size_t activeCellCount = std::count_if(
        cells.begin(),
        cells.end(),
        [](const CellData& cell) { return !cell.isDiscarded; });

    parallelFor(allCells, [&](uint idx) {
        if (cells[idx].isDiscarded) {
            return;
        }

        cells[idx] = shootRayOnCell(
            cells[idx],
            m,
            scene,
            planePoint,
            direction,
            MAX_DISTANCE,
            RAY_EPS);
    });

    uint debugStepIndex = 1;
    if (debug) {
        saveMoldCheckStepMesh( // Step 1
            cells,
            direction,
            debugResultsSubdir,
            debugStepIndex);
    }

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
                  << rawHitCount << "/" << activeCellCount << "\n";
        std::cout << "Beginning Clamping phase...\n";
        std::cout.flush();
    }

    parallelFor(allCells, [&](uint idx) {
        if (cells[idx].isDiscarded) {
            return;
        }

        computeClampedCell(
            idx,
            cells,
            hitCellIds,
            planePoint,
            direction,
            CONE_COS_THRESHOLD,
            EPS);
    });

    if (debug) {
        saveMoldCheckStepMesh( // Step 2
            cells,
            direction,
            debugResultsSubdir,
            debugStepIndex);
    }

    const double REDUCE_POINTS_ANGLE_THRESHOLD_DEGREES = 10.0;
    
    cells = reducePoints(
        cells,
        grid,
        direction,
        draftAngleDegrees,
        EPS,
        REDUCE_POINTS_ANGLE_THRESHOLD_DEGREES,
        MAX_DISTANCE,
        debug,
        debugResultsSubdir,
        &debugStepIndex);

    if (hasEnclosedWhiteHole(
            cells,
            grid,
            direction,
            debug,
            debugResultsSubdir)) {
        std::cout << "Direction discarded: the reduced cell grid contains "
                     "an enclosed hole of white cells.\n";
        std::cout.flush();
        return {};
    }

    double totalAreaHit = 0.0;
    double clampedAreaHit = 0.0;
    double hiddenAreaHit = 0.0;
    uint reducedHitCount = 0;

    for (uint i = 0; i < cells.size(); ++i) {
        if (cells[i].hasHit) {
            ++reducedHitCount;
            totalAreaHit += cellArea;
        }

		// if (cells[i].hasHit &&
		// 	cells[i].hasClampedHit &&
		// 	std::abs(cells[i].clampedDistance - cells[i].distance) > EPS) {
		if (cells[i].hasHit && cells[i].hasClampedHit) {
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
        (activeCellCount > 0) ?
            static_cast<double>(reducedHitCount) / activeCellCount :
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
                  << reducedHitCount << "/" << activeCellCount << "\n";
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
                magentaAngleDegrees,
                magentaCellInterval,
                MAX_BIHARMONIC_ITERATIONS,
                EPS,
                debugResultsSubdir,
                MAX_DISTANCE,
                &debugStepIndex);

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
				// cells[i].cellCenter + direction * cells[i].clampedDistance,
				cells[i].cellCenter + direction * cells[i].distance,
                Color::Blue);
        }

        PolyMesh depthPointsMesh;
        depthPointsMesh.enablePerVertexColor();

        for (uint i = 0; i < depthCells.size(); ++i) {
            if (depthCells[i].isDiscarded) {
                continue;
            }

            const Point3d depthPoint =
                depthCells[i].cellCenter +
                direction * depthCells[i].distance;

            addColoredPoint(
                depthPointsMesh,
                depthPoint,
                moldCheckCellDebugColor(depthCells[i]));
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
#endif
