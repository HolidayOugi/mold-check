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
#include <vclib/io.h>
#include <vclib/meshes.h>
#include <vclib/algorithms/core/fibonacci.h>



#include <chrono>
#include <cmath>
#include <filesystem>
#include <limits>
#include <numeric>
#include <unordered_set>
#include <vector>

#include "include/struct.h"
#include "include/helper.h"
#include "include/functions.h"
#include "include/debug_output.h"

struct MoldCheckMetrics
{
	double score = -std::numeric_limits<double>::infinity();
	double componentRatio = 0.0;
	vcl::uint largestComponentSize = 0;
	double percentClamped = 0.0;
	double percentHidden = 0.0;
};

static double moldQualityScore(
	double componentRatio,
	double percentClamped,
	double percentHidden)
{
	return
		0.50 * componentRatio +
		0.30 * (1 - (percentClamped / 100.0)) +
		0.20 * (1 - (percentHidden / 100.0));
}

MoldCheckMetrics moldCheck(
	vcl::PolyMesh              m,
	const std::vector<double>& gridCellSideLengths,
	bool                       debug,
	vcl::Point3d 			   direction,
	const double 			   coneAngleDegrees,
	const double 			   marginFactor,
	const std::string&         debugResultsSubdir = "",
	const vcl::PolyMesh*       mold = nullptr,
	vcl::TriMesh*              outRemainingMoldMesh = nullptr)
{
	using namespace vcl;

	const double CONE_COS_THRESHOLD = std::cos(coneAngleDegrees * M_PI / 180.0);
    
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
	const Planed  plane(planePoint, direction);
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

    const double cellDu   = grid.sideU;
	const double cellDv   = grid.sideV;
	const double cellArea = cellDu * cellDv;
	std::vector<uint> allCells(grid.rows * grid.cols);
	std::iota(allCells.begin(), allCells.end(), 0);
	std::vector<CellData> cells(allCells.size());

	parallelFor(allCells, [&](uint idx) {
		const CellData cell = makeCellGeometry(idx, grid, planePoint, u, v);
		cells[idx] = shootRayOnCell(cell, m, scene, planePoint, direction, MAX_DISTANCE, RAY_EPS);
	});

    std::vector<CellData> clampedCells = cells;

        
	if (debug) {
            std::cout << "Ray casting complete. Hit points: ";
            uint hitCount = 0;
			for (uint i = 0; i < cells.size(); ++i) {
				if (cells[i].distance < MAX_DISTANCE) {
                    ++hitCount;
                }
            }
            std::cout << hitCount << "/" << allCells.size() << " cells \n";
			std::cout << "Beginning Clamping phase...\n";
			std::cout.flush();
	}


	parallelFor(allCells, [&](uint idx) {
		clampedCells[idx] = computeClampedCell(idx, cells, planePoint, direction, CONE_COS_THRESHOLD, EPS);
	});

	double totalAreaHit = 0.0;
	double clampedAreaHit = 0.0;
	double hiddenAreaHit = 0.0;

	for (uint i = 0; i < cells.size(); ++i) {
		if (!cells[i].hasHit) {
			continue;
		}

		totalAreaHit += cellArea;

		if (clampedCells[i].distance != cells[i].distance) {
			clampedAreaHit += cellArea;
		}

		if (cells[i].hitPoints.size() > 2) {
			hiddenAreaHit += cellArea;
		}
	}

	const double percentClamped =
		(totalAreaHit > 0.0) ? (clampedAreaHit / totalAreaHit) * 100.0 : 0.0;
	const double percentHidden =
		(totalAreaHit > 0.0) ? (hiddenAreaHit / totalAreaHit) * 100.0 : 0.0;

	ConnectedComponentData largestComponent =
		largestConnectedComponent(
			cells, clampedCells, grid, EPS);

	const double componentRatio =
		(totalAreaHit > 0.0) ? (largestComponent.area / totalAreaHit) : 0.0;

	const MoldCheckMetrics metrics{
		moldQualityScore(componentRatio, percentClamped, percentHidden),
		componentRatio,
		static_cast<uint>(largestComponent.indices.size()),
		percentClamped,
		percentHidden};

	std::vector<CellData> moldClampedCells = clampedCells;

	if (mold != nullptr) {
		embree::Scene moldScene(*mold);
		for (uint i = 0; i < moldClampedCells.size(); ++i) {
			if (!moldClampedCells[i].hasHit) continue;

			const Point3d rayOrigin =
				moldClampedCells[i].hitPoints[0] - direction * RAY_EPS;

			const auto rayHits =
				moldScene.facesIntersectedByRay(
					rayOrigin,
					direction,
					RAY_EPS);

			if (rayHits.empty()) {
				moldClampedCells[i].hasHit = false;
				continue;
			};

			const auto [faceId, baryCoords, triId, hitT] =
				rayHits.back();

			moldClampedCells[i].hitPoints = {
				computeHitPoint(
					*mold,
					faceId,
					triId,
					baryCoords,
					moldClampedCells[i].hitPoints[0])};
			moldClampedCells[i].distance = hitT;
		}
	}

	if (outRemainingMoldMesh != nullptr) {
		*outRemainingMoldMesh =
			createRemainingMold(
				cells,
				clampedCells,
				moldClampedCells,
				direction);
	}

	

	if (debug) {
		std::cout << "Validating clamped cells...\n";
		std::cout.flush();

		validateClampedCells(clampedCells, allCells, direction, CONE_COS_THRESHOLD, EPS);
		
		debugOutput(
			cells,
			clampedCells,
			moldClampedCells,
			largestComponent,
			grid,
			planePoint,
			u,
			v,
			direction,
			EPS,
			debugResultsSubdir,
			totalAreaHit,
			clampedAreaHit,
			percentClamped,
			hiddenAreaHit,
			percentHidden,
			componentRatio,
			metrics.score);
    }
        
    return metrics;
}

int main()
{
    using namespace vcl;

	const auto startTime = std::chrono::steady_clock::now();

	const uint NUM_PLANES = 10;

	std::vector<Point3d> fibNormals = sphericalFibonacciPointSet<Point3d>(NUM_PLANES);


    PolyMesh m = loadMesh<PolyMesh>(MESHES_PATH "/bimba_enlarged.ply");


    std::vector<double> gridCellSideLengths = {0.3, 0.3};

	const double coneAngleDegrees = 5.0;

	const double marginFactor = 0.1;

	PolyMesh mold = squareMold(m, 0.05);
	const std::filesystem::path externalResultsPath =
		RESULTS_PATH;
	std::filesystem::create_directories(externalResultsPath);
	saveMesh(
		mold,
		(externalResultsPath / "mold.ply").string());
	
	double moldVolume = squareMoldVolume(m, marginFactor);
	std::cout << "Mold volume: " << moldVolume << "\n";

	MoldCheckMetrics result;
	MoldCheckMetrics bestResult;
	MoldCheckMetrics worstResult;
	worstResult.score = std::numeric_limits<double>::infinity();
	int bestDirectionIndex = 0;
	int worstDirectionIndex = 0;

	for (const auto& direction : fibNormals) {
		std::cout << "Processing direction: " << direction << "\n";
		result = moldCheck(m, gridCellSideLengths, false, direction, coneAngleDegrees, marginFactor);
		std::cout << "Score: " << result.score
					  << " (component ratio: " << result.componentRatio
					  << ", largest component cells: " << result.largestComponentSize
					  << ", clamped: " << result.percentClamped << "%"
					  << ", hidden: " << result.percentHidden << "%)\n";
		if (result.score > bestResult.score) {
			bestResult = result;
			bestDirectionIndex = &direction - &fibNormals[0];
			std::cout << ("New best direction found! Index: " + std::to_string(bestDirectionIndex) + ", Score: " + std::to_string(bestResult.score) + "\n");
		}
		if (result.score < worstResult.score) {
			worstResult = result;
			worstDirectionIndex = &direction - &fibNormals[0];
			std::cout << ("New worst direction found! Index: " + std::to_string(worstDirectionIndex) + ", Score: " + std::to_string(worstResult.score) + "\n");
		}
	}

	TriMesh bestRemainingMoldMesh;
	result = moldCheck(
		m,
		gridCellSideLengths,
		true,
		fibNormals[bestDirectionIndex],
		coneAngleDegrees,
		marginFactor,
		"best",
		&mold,
		&bestRemainingMoldMesh);
	const std::filesystem::path bestRemainingMoldPath =
		externalResultsPath /
		"best" / "mold" / "888_mold_check_remaining_mold.ply";
	std::filesystem::create_directories(bestRemainingMoldPath.parent_path());
	saveMesh(bestRemainingMoldMesh, bestRemainingMoldPath.string());

	TriMesh worstRemainingMoldMesh;
	result = moldCheck(
		m,
		gridCellSideLengths,
		true,
		fibNormals[worstDirectionIndex],
		coneAngleDegrees,
		marginFactor,
		"worst",
		&mold,
		&worstRemainingMoldMesh);
	const std::filesystem::path worstRemainingMoldPath =
		externalResultsPath /
		"worst" / "mold" / "888_mold_check_remaining_mold.ply";
	std::filesystem::create_directories(worstRemainingMoldPath.parent_path());
	saveMesh(worstRemainingMoldMesh, worstRemainingMoldPath.string());

    const auto endTime = std::chrono::steady_clock::now();
    const auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "moldCheck execution time: " << elapsedMs.count() << " ms\n";
    std::cout.flush();
    return 0;
}
