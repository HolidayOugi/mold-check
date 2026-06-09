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



#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <limits>
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
	double percentClamped = 0.0;
	double percentHidden = 0.0;
};

static double moldQualityScore(
	double hitRatio,
	double compactness,
	double percentClamped,
	double percentHidden)
{
	return
		0.35 * hitRatio +
		0.25 * compactness +
		0.20 * (1 - (percentHidden / 100.0)) +
		0.20 * (1 - (percentClamped / 100.0));
}

MoldCheckMetrics moldCheck(
	vcl::PolyMesh              m,
	const std::vector<double>& gridCellSideLengths,
	bool                       debug,
	vcl::Point3d 			   direction,
	const double 			   coneAngleDegrees,
	const double 			   marginFactor,
	const std::string&         debugResultsSubdir = "")
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

	uint rawHitCount = 0;
        
	if (debug) {
			for (uint i = 0; i < cells.size(); ++i) {
				if (cells[i].hasHit) {
                    ++rawHitCount;
                }
			}
			std::cout << "Ray casting complete. Hit cells: "
					  << rawHitCount << "/" << allCells.size() << "\n";
			std::cout << "Beginning Clamping phase...\n";
			std::cout.flush();
	}


	parallelFor(allCells, [&](uint idx) {
		computeClampedCell(idx, cells, planePoint, direction, CONE_COS_THRESHOLD, EPS);
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
		REDUCE_POINTS_DISTANCE_THRESHOLD);


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
		(totalAreaHit > 0.0) ? (clampedAreaHit / totalAreaHit) * 100.0 : 0.0;
	const double percentHidden =
		(totalAreaHit > 0.0) ? (hiddenAreaHit / totalAreaHit) * 100.0 : 0.0;

	const HitCellShapeData hitShape = hitCellShape(cells, grid);

	const double hitRatio =
		(cells.size() > 0) ?
			static_cast<double>(reducedHitCount) / cells.size() :
			0.0;

	const MoldCheckMetrics metrics{
		moldQualityScore(
			hitRatio,
			hitShape.compactness,
			percentClamped,
			percentHidden),
		hitRatio,
		hitShape.compactness,
		reducedHitCount,
		percentClamped,
		percentHidden};

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

		PolyMesh violatingPointsMesh = validateClampedCells(depthCells, allCells, direction, CONE_COS_THRESHOLD, EPS);
		
		PolyMesh hitPointsMesh;
		hitPointsMesh.enablePerVertexColor();
		for (uint i = 0; i < cells.size(); ++i) {
			if (cells[i].distance == MAX_DISTANCE) continue;
			addColoredPoint(
				hitPointsMesh,
				cells[i].cellCenter + direction * cells[i].distance,
				Color::Yellow);
		}

		PolyMesh hitPointsafterReductionMesh;
		hitPointsafterReductionMesh.enablePerVertexColor();
		for (uint i = 0; i < cells.size(); ++i) {
			if (!cells[i].hasHit) continue;
			addColoredPoint(
				hitPointsafterReductionMesh,
				cells[i].cellCenter + direction * cells[i].distance,
				Color::Blue);
		}
		
		PolyMesh clampedPointsMesh;
		clampedPointsMesh.enablePerVertexColor();
		for (uint i = 0; i < cells.size(); ++i) {
			if (!cells[i].hasHit) continue;
			addColoredPoint(
				clampedPointsMesh,
				cells[i].cellCenter + direction * cells[i].clampedDistance,
				Color::Blue);
		}

		PolyMesh depthPointsMesh;
		depthPointsMesh.enablePerVertexColor();
		for (uint i = 0; i < depthCells.size(); ++i) {
			const Point3d depthPoint =
				depthCells[i].cellCenter + direction * depthCells[i].distance;
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

		const TriMesh moldSurfaceMesh = createMoldSurface(depthCells, grid, direction);

		
		const std::filesystem::path debugOutputDir =
			std::filesystem::path(RESULTS_PATH) /
			debugResultsSubdir;

		std::filesystem::create_directories(debugOutputDir);

		/*
		for (const auto& entry :
			 std::filesystem::directory_iterator(debugOutputDir)) {
			if (entry.is_regular_file() && entry.path().extension() == ".ply") {
				std::filesystem::remove(entry.path());
			}
		}

		*/

		const std::string base =
			(debugOutputDir / "mold_check").string();
		saveMesh(hitPointsMesh, base + "_hit_points.ply");
		saveMesh(hitPointsafterReductionMesh, base + "_hit_points_after_reduction.ply");
		saveMesh(clampedPointsMesh, base + "_lipschitz_points.ply");
		saveMesh(depthPointsMesh, base + "_mold_points.ply");
		saveMesh(planeMesh, base + "_plane_tested.ply");
		saveMesh(moldSurfaceMesh, base + "_mold_surface.ply");
		saveMesh(violatingPointsMesh, base + "_non-lipschitz_points.ply");
		
		std::cout << "Clamped points: " << clampedPointsMesh.vertexCount() << "\n";
		std::cout << "Depth points: " << depthPointsMesh.vertexCount() << "\n";
		std::cout << "Mold surface median points: " << moldSurfaceMesh.vertexCount() << "\n";
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
		std::cout << "percentHidden: "
				  << percentHidden << "\n";
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
				<< " - " << base << "_non-lipschitz_points.ply\n";

		std::cout << "=== moldCheck completed successfully ===\n";
		std::cout.flush();
    }
        
    return metrics;
}

int main()
{
    using namespace vcl;

	const auto startTime = std::chrono::steady_clock::now();

	const uint NUM_PLANES = 100;

	std::vector<Point3d> fibNormals = sphericalFibonacciPointSet<Point3d>(NUM_PLANES);


    PolyMesh m = loadMesh<PolyMesh>(MESHES_PATH "/bimba_enlarged.ply");


    std::vector<double> gridCellSideLengths = {0.4, 0.4};

	const double coneAngleDegrees = 5.0;

	const double marginFactor = 0.05;

	const std::filesystem::path externalResultsPath =
		RESULTS_PATH;
	std::filesystem::create_directories(externalResultsPath);

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
		std::cout << "Processing direction: " << direction << "\n";
		result = moldCheck(m, gridCellSideLengths, false, direction, coneAngleDegrees, marginFactor);
		scoredDirections.push_back({result.score, directionIndex});
		std::cout << "Score: " << result.score
					  << " (hit ratio: " << result.hitRatio
					  << ", compactness: " << result.compactness
					  << ", hits: " << result.hitCount
					  << ", clamped: " << result.percentClamped << "%"
					  << ", hidden: " << result.percentHidden << "%)\n";
		if (result.score > bestResult.score) {
			bestResult = result;
			bestDirectionIndex = directionIndex;
			std::cout << ("New best direction found! Index: " + std::to_string(bestDirectionIndex) + ", Score: " + std::to_string(bestResult.score) + "\n");
		}
		if (result.score < worstResult.score) {
			worstResult = result;
			worstDirectionIndex = directionIndex;
			std::cout << ("New worst direction found! Index: " + std::to_string(worstDirectionIndex) + ", Score: " + std::to_string(worstResult.score) + "\n");
		}
	}

	std::sort(
		scoredDirections.begin(),
		scoredDirections.end(),
		[](const auto& a, const auto& b) {
			return a.first < b.first;
		});

	result = moldCheck(
		m,
		gridCellSideLengths,
		true,
		fibNormals[bestDirectionIndex],
		coneAngleDegrees,
		marginFactor,
		"best");

	result = moldCheck(
		m,
		gridCellSideLengths,
		true,
		fibNormals[worstDirectionIndex],
		coneAngleDegrees,
		marginFactor,
		"worst");

	const int medianDebugCount =
		std::min<int>(15, static_cast<int>(scoredDirections.size()));
	const int medianCenter =
		static_cast<int>(scoredDirections.size() / 2);
	const int medianStart =
		std::max(
			0,
			std::min(
				medianCenter - medianDebugCount / 2,
				static_cast<int>(scoredDirections.size()) - medianDebugCount));

	for (int i = 0; i < medianDebugCount; ++i) {
		const auto& [medianScore, medianDirectionIndex] =
			scoredDirections[medianStart + i];
		const std::string medianDir =
			"median " + std::to_string(i + 1);

		std::cout << "Processing median direction " << (i + 1)
				  << "/" << medianDebugCount
				  << ". Index: " << medianDirectionIndex
				  << ", Score: " << medianScore
				  << ", Output: " << medianDir << "\n";

		result = moldCheck(
			m,
			gridCellSideLengths,
			true,
			fibNormals[medianDirectionIndex],
			coneAngleDegrees,
			marginFactor,
			medianDir);
	}

    const auto endTime = std::chrono::steady_clock::now();
    const auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "moldCheck execution time: " << elapsedMs.count() << " ms\n";
    std::cout.flush();
    return 0;
}
