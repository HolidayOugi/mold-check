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

#include "include/mold_check.h"

#include <vclib/algorithms/core/fibonacci.h>
#include <vclib/algorithms/mesh/update/bounding_box.h>
#include <vclib/igl/booleans.h>
#include <vclib/io.h>
#include <vclib/meshes.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

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

    const std::filesystem::path resultsPath = RESULTS_PATH;
    try {
        std::filesystem::create_directories(resultsPath);
        for (const auto& entry : std::filesystem::directory_iterator(resultsPath)) {
            std::filesystem::remove_all(entry.path());
        }
    }
    catch (const std::filesystem::filesystem_error& error) {
        std::cerr << "Error: unable to clear output directory '"
                  << resultsPath.string() << "': " << error.what() << "\n";
    }

    //bigger than boxMesh for it to completely encapsulate the mold mesh
	const double marginFactor = 0.4;
    const uint NUM_PLANES = 100;

    std::vector<Point3d> fibNormals =
        sphericalFibonacciPointSet<Point3d>(NUM_PLANES);

    PolyMesh m =
        loadMesh<PolyMesh>(MESHES_PATH "/bimba_enlarged.ply");

    const PolyMesh moldMesh =
        makeContainingBoxMesh(m, 0.2);

    saveMesh(moldMesh, RESULTS_PATH "/mold-bool.ply");

    std::vector<double> gridCellSideLengths = {0.4, 0.4};

    const double coneAngleDegrees = 5.0;
    const double magentaAngleDegrees = 45.0;
    const vcl::uint magentaCellInterval = 20;
    const double draftAngleDegrees = 20.0;

    MoldCheckMetrics result;
    MoldCheckMetrics bestResult;
    MoldCheckMetrics worstResult;

    worstResult.score = std::numeric_limits<double>::infinity();

    int bestDirectionIndex = 0;
    int worstDirectionIndex = 0;

    std::vector<std::pair<double, int>> scoredDirections;

    //Temporarily test only direction 0 instead of all 100 directions.
    const std::vector<int> directionIndicesToTest = {0};

    //for (const int directionIndex : directionIndicesToTest) {
    for (uint directionIndex = 0; directionIndex < fibNormals.size(); ++directionIndex) {
        const auto& direction = fibNormals[directionIndex];

        std::cout << "Processing direction: "
                  << direction << "\n";

        result =
            moldCheck(
                m,
                gridCellSideLengths,
                false,
                direction,
                coneAngleDegrees,
                magentaAngleDegrees,
                magentaCellInterval,
                draftAngleDegrees,
                marginFactor);

        if (!std::isfinite(result.score)) {
            continue;
        }

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
                magentaAngleDegrees,
                magentaCellInterval,
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
            magentaAngleDegrees,
            magentaCellInterval,
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
            magentaAngleDegrees,
            magentaCellInterval,
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
                magentaAngleDegrees,
                magentaCellInterval,
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
