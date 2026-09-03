#ifndef VCL_TEST_EXTERNAL_888_MOLD_CHECK_FUNCTIONS_H
#define VCL_TEST_EXTERNAL_888_MOLD_CHECK_FUNCTIONS_H

#include "helper.h"
#include "biharmonic.h"
#include "struct.h"
#include "debug_output.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <unordered_set>
#include <vector>

#include <vclib/algorithms/mesh/update/bounding_box.h>
#include <vclib/embree/scene.h>
#include <vclib/meshes.h>



static std::vector<CellData> fixDepthCellConeViolations(
	std::vector<CellData> depthCells,
	const vcl::Point3d& direction,
	double coneCosThreshold,
	float eps)
{
	using namespace vcl;

	const std::vector<CellData> originalDepthCells = depthCells;
	std::vector<uint> allCells(depthCells.size());
	std::iota(allCells.begin(), allCells.end(), 0);

	parallelFor(allCells, [&](uint i) {
		if (originalDepthCells[i].isDiscarded ||
			originalDepthCells[i].hitPoints.empty()) {
			return;
		}

		const Point3d original = originalDepthCells[i].hitPoints[0];

		double requiredT = 0.0;
		bool hasViolation = false;

		for (uint j = 0; j < originalDepthCells.size(); ++j) {
			if (i == j ||
				originalDepthCells[j].isDiscarded ||
				originalDepthCells[j].hitPoints.empty()) {
				continue;
			}

			if (!isWithinPlaneAngle(
					original,
					originalDepthCells[j].hitPoints[0],
					direction,
					coneCosThreshold,
					eps)) {
				continue;
			}

			hasViolation = true;

			const double t = coneBoundaryStep(
				original,
				originalDepthCells[j].hitPoints[0],
				direction,
				coneCosThreshold,
				eps);

			requiredT = std::max(requiredT, t);
		}

		if (!hasViolation) {
			return;
		}

		const Point3d fixedPoint = original - direction * requiredT;

		depthCells[i].hitPoints[0] = fixedPoint;
		depthCells[i].distance =
			(fixedPoint - depthCells[i].cellCenter).dot(direction);
		// depthCells[i].clampedDistance = depthCells[i].distance;
		depthCells[i].hasClampedHit = true;
	});

	return depthCells;
}

static std::vector<CellData> makeDepthCells(
	const std::vector<CellData>& cells,
	const vcl::Point3d& direction,
	const GridChoice& grid,
	double coneCosThreshold,
	float eps,
	const std::string&         debugResultsSubdir = "",
	double maxDistance = std::numeric_limits<double>::infinity(),
	vcl::uint* debugStepIndex = nullptr)
{
	using namespace vcl;

	if (cells.size() != grid.rows * grid.cols) {
		return {};
	}


	if (std::none_of(cells.begin(), cells.end(), [](const CellData& cell) {
			return cell.hasHit;
		})) {
		return cells;
	}
	std::vector<CellData> surfaceCells = cells;
	std::vector<CellData> depthCells = surfaceCells;
	updateDepthCellInsideFlags(surfaceCells, depthCells, eps);

	depthCells =
		biharmonicFillWhiteCells(
			surfaceCells,
			depthCells,
			grid,
			direction,
			eps);

	if (debugStepIndex != nullptr) {
		saveMoldCheckStepMesh( // Step 10
			depthCells,
			direction,
			debugResultsSubdir,
			*debugStepIndex);
	}

	PolyMesh depthPointsMesh;
	depthPointsMesh.enablePerVertexColor();
	for (uint i = 0; i < depthCells.size(); ++i) {
		if (depthCells[i].isDiscarded) {
			continue;
		}

		const Point3d depthPoint =
			depthCells[i].cellCenter + direction * depthCells[i].distance;
		addColoredPoint(
			depthPointsMesh,
			depthPoint,
			moldCheckCellDebugColor(depthCells[i]));
	}

	const std::filesystem::path debugOutputDir =
			std::filesystem::path(RESULTS_PATH) /
			debugResultsSubdir;

	std::filesystem::create_directories(debugOutputDir);
	
	const std::string base =
			(debugOutputDir / "mold_check").string();
	saveMesh(depthPointsMesh, base + "_original_mold_points.ply");

	depthCells = biharmonicFillHitCells(
		surfaceCells,
		depthCells,
		grid,
		direction,
		eps,
		0,
		maxDistance);

	if (debugStepIndex != nullptr) {
		saveMoldCheckStepMesh( // Step 11
			depthCells,
			direction,
			debugResultsSubdir,
			*debugStepIndex);
	}

	depthCells =
		biharmonicFillWhiteCells(
			surfaceCells,
			depthCells,
			grid,
			direction,
			eps,
			maxDistance);
	
	if (debugStepIndex != nullptr) {
		saveMoldCheckStepMesh( // Step 12
			depthCells,
			direction,
			debugResultsSubdir,
			*debugStepIndex);
	}

	// Remember which white cells entered this pass as forward-capped cyan.
	std::vector<unsigned char> cyanCells(depthCells.size(), false);
	for (uint idx = 0; idx < surfaceCells.size(); ++idx) {
		cyanCells[idx] =
			biharmonicIsWhiteForwardCapCandidate(
				surfaceCells,
				depthCells,
				idx) &&
			depthCells[idx].isMovedForward;
	}

	depthCells = biharmonicFillHitCells(
		surfaceCells,
		depthCells,
		grid,
		direction,
		eps,
		10,
		maxDistance,
		&cyanCells);
	
	if (debugStepIndex != nullptr) {
		saveMoldCheckStepMesh( // Step 13
			depthCells,
			direction,
			debugResultsSubdir,
			*debugStepIndex);
	}

	depthCells = fixDepthCellConeViolations(depthCells, direction, coneCosThreshold, eps);
	updateDepthCellInsideFlags(surfaceCells, depthCells, eps);

	if (debugStepIndex != nullptr) {
		saveMoldCheckStepMesh( // Step 14
			depthCells,
			direction,
			debugResultsSubdir,
			*debugStepIndex);
	}


	return depthCells;
}

#endif
