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



static std::tuple<vcl::Point3d, vcl::Point3d> makePlane(
	const vcl::PolyMesh& m,
	const vcl::Planed&   plane,
	const vcl::Point3d&  planePoint,
	const vcl::Point3d&  direction,
	double               margin,
	double               eps,
	GridChoice&          grid)
{
	using namespace vcl;

	Point3d u;
	Point3d v;

	direction.orthoBase(u, v);

	if (u.norm() <= eps || v.norm() <= eps) {
		grid.minU = 0.0;
		grid.minV = 0.0;
		grid.maxU = 0.0;
		grid.maxV = 0.0;
		return {u, v};
	}

	u.normalize();
	v.normalize();

	grid.minU = std::numeric_limits<double>::infinity();
	grid.minV = std::numeric_limits<double>::infinity();
	grid.maxU = -std::numeric_limits<double>::infinity();
	grid.maxV = -std::numeric_limits<double>::infinity();

	for (const auto& vert : m.vertices()) {
		const Point3d projected = plane.projectPoint(vert.position());
		const Point3d rel       = projected - planePoint;

		const double pu = rel.dot(u);
		const double pv = rel.dot(v);

		grid.minU = std::min(grid.minU, pu - margin);
		grid.minV = std::min(grid.minV, pv - margin);
		grid.maxU = std::max(grid.maxU, pu + margin);
		grid.maxV = std::max(grid.maxV, pv + margin);
	}

	return {u, v};
}

static void makeGrid(
	GridChoice&                grid,
	const std::vector<double>& gridCellSideLengths)
{
	using namespace vcl;

	const double lenU = grid.maxU - grid.minU;
	const double lenV = grid.maxV - grid.minV;

	const double sideU =
		(gridCellSideLengths.size() >= 1) ? gridCellSideLengths[0] : lenU;

	const double sideV =
		(gridCellSideLengths.size() >= 2) ? gridCellSideLengths[1] : sideU;

	if (sideU <= 0.0 || sideV <= 0.0) {
		grid.rows = 1;
		grid.cols = 1;
		grid.sideU = lenU;
		grid.sideV = lenV;
		return;
	}

	grid.cols = static_cast<uint>(std::max(1.0, std::ceil(lenU / sideU)));
	grid.rows = static_cast<uint>(std::max(1.0, std::ceil(lenV / sideV)));

	grid.sideU = sideU;
	grid.sideV = sideV;

	grid.maxU = grid.minU + grid.cols * grid.sideU;
	grid.maxV = grid.minV + grid.rows * grid.sideV;
}

static CellData shootRayOnCell(
	const CellData& cell,
	const vcl::PolyMesh& m,
	const vcl::embree::Scene& scene,
	const vcl::Point3d& planePoint,
	const vcl::Point3d& direction,
	double maxDistance,
	float eps)
{
	using namespace vcl;

	const Point3d rayOrigin =
		cell.cellCenter + direction * (-eps);

	const Point3d invalidPoint =
		cell.cellCenter + direction * maxDistance;

	//noticeably faster to use firstFaceIntersectedbyRay and then recast multiple
	//rays only on non-empty cells than to use facesIntersectedByRay directly for all cells
	auto [faceId, baryCoords, triId, hitT] =
		scene.firstFaceIntersectedByRay(rayOrigin, direction);

	if (faceId != UINT_NULL) {
		//redoing the first hit might seem redudant but it is actually faster than computing the hit point three times.
		//facesIntersectedbyRay has -eps built in
		auto rayHits = scene.facesIntersectedByRay(rayOrigin, direction, eps);
		std::unordered_set<uint> seenFaceIds;
		decltype(rayHits) uniqueRayHits;
		uniqueRayHits.reserve(rayHits.size());

		for (const auto& rayHit : rayHits) {
			const uint hitFaceId = std::get<0>(rayHit);
			if (seenFaceIds.insert(hitFaceId).second) {
				uniqueRayHits.push_back(rayHit);
			}
		}

		rayHits = std::move(uniqueRayHits);

		//fallback for possible missed hit due to numerical issues in firstFaceIntersectedByRay and silent crash
		//possible bug to fix
		if (rayHits.empty()) {
			CellData result = cell;
			result.distance = hitT;
			result.hitPoints = {
				computeHitPoint(m, faceId, triId, baryCoords, invalidPoint)};
			result.hasHit = result.hitPoints[0] != invalidPoint;
			result.hasClampedHit = false;
			result.clampedDistance = result.distance;

			return result;
		}
		
		auto [faceId, baryCoords, triId, hitT] = rayHits.front(); 
		Point3d hitPoint = computeHitPoint(m, faceId, triId, baryCoords, invalidPoint);

		if (hitPoint != invalidPoint) {
			CellData result = cell;
			result.distance = hitT;
			result.hitPoints.clear();
			result.hitPoints.reserve(rayHits.size());
			for (const auto& rayHit : rayHits) {
				result.hitPoints.push_back(
					computeHitPoint(
						m,
						std::get<0>(rayHit),
						std::get<2>(rayHit),
						std::get<1>(rayHit),
						invalidPoint));
			}
			result.hasHit = true;
			result.hasClampedHit = false;
			result.clampedDistance = result.distance;

			return result;
		}
	
	}

	CellData result = cell;
	result.distance = maxDistance;
	result.hitPoints = {invalidPoint};
	result.hasHit = false;
	result.hasClampedHit = false;
	result.clampedDistance = maxDistance;

	return result;
}

static CellData makeCellGeometry(
	vcl::uint idx,
	const GridChoice& grid,
	const vcl::Point3d& planePoint,
	const vcl::Point3d& u,
	const vcl::Point3d& v)
{
	using namespace vcl;

	const uint j = idx / grid.cols;
	const uint i = idx % grid.cols;

	const double u0 = grid.minU + i * grid.sideU;
	const double u1 = u0 + grid.sideU;

	const double v0 = grid.minV + j * grid.sideV;
	const double v1 = v0 + grid.sideV;

	CellData cell;

	cell.cellCorners = {
		planePoint + u * u0 + v * v0,
		planePoint + u * u1 + v * v0,
		planePoint + u * u1 + v * v1,
		planePoint + u * u0 + v * v1};

	const double centerU = grid.minU + (i + 0.5) * grid.sideU;
	const double centerV = grid.minV + (j + 0.5) * grid.sideV;

	cell.cellCenter = planePoint + u * centerU + v * centerV;

	cell.distance = 0.0;
	cell.hitPoints = {cell.cellCenter};
	cell.hasHit = false;
	cell.hasClampedHit = false;
	cell.clampedDistance = 0.0;

	return cell;
}

static void computeClampedCell(
	vcl::uint i,
	std::vector<CellData>& cells,
	const std::vector<vcl::uint>& hitCellIds,
	const vcl::Point3d& planePoint,
	const vcl::Point3d& direction,
	double coneCosThreshold,
	float eps)
{
	using namespace vcl;

	const CellData baseCell = cells[i];

	if (!baseCell.hasHit) {
		cells[i].hasClampedHit = false;
		cells[i].clampedDistance = baseCell.distance;
		return;
	}

	const Point3d original = baseCell.hitPoints[0];

	double requiredT = 0.0;
	bool anyCone = false;

	for (uint j : hitCellIds) {
		if (i == j) {
			continue;
		}

		if (!isWithinPlaneAngle(
				original,
				cells[j].hitPoints[0],
				direction,
				coneCosThreshold,
				eps)) {
			continue;
		}

		anyCone = true;

		const double t = coneBoundaryStep(
			original,
			cells[j].hitPoints[0],
			direction,
			coneCosThreshold,
			eps);

		requiredT = std::max(requiredT, t);
	}

	if (!anyCone) {
		cells[i].hasClampedHit = false;
		cells[i].clampedDistance = baseCell.distance;
		return;
	}

	const Point3d currentPoint =
		original - direction * requiredT;

	const double distanceToPlane =
		std::abs((currentPoint - planePoint).dot(direction));

	cells[i].hasClampedHit = true;
	cells[i].clampedDistance = distanceToPlane;
}

static void removeDraftAngleBoundaryPoints(
	std::vector<CellData>& cells,
	const GridChoice& grid,
	const vcl::Point3d& direction,
	double draftAngleDegrees,
	float eps,
	double maxDistance)
{
	using namespace vcl;

	if (cells.size() != grid.rows * grid.cols ||
		draftAngleDegrees <= 0.0 || direction.norm() == 0.0) {
		return;
	}

	Point3d normalizedDirection = direction;
	normalizedDirection.normalize();
	const double coneCosThreshold =
		std::cos(draftAngleDegrees * M_PI / 180.0);

	std::vector<uint> allCells(cells.size());
	std::iota(allCells.begin(), allCells.end(), 0);
	std::vector<char> removeCell(cells.size(), false);
	std::atomic<bool> removedAny = false;

	struct DraftAnglePoint
	{
		uint idx = 0;
		Point3d point;
	};

	do {
		removedAny.store(false, std::memory_order_relaxed);
		std::fill(removeCell.begin(), removeCell.end(), false);

		//precompute the points to compare against for efficiency
		std::vector<DraftAnglePoint> comparisonPoints;
		comparisonPoints.reserve(cells.size());
		for (uint otherIdx = 0; otherIdx < cells.size(); ++otherIdx) {
			const CellData& other = cells[otherIdx];
			if ((!other.hasHit && !other.hasClampedHit) ||
				other.hitPoints.empty()) {
				continue;
			}

			comparisonPoints.push_back(
				{
					otherIdx,
					other.hasClampedHit ?
						other.cellCenter +
							normalizedDirection * other.clampedDistance :
						other.hitPoints[0]
				});
		}

		parallelFor(allCells, [&](uint idx) {
			if (!cells[idx].hasHit || cells[idx].hasClampedHit ||
				cells[idx].hitPoints.empty()) {
				return;
			}

			bool isBoundary = false;
			forEachSquareNeighbor(idx, grid, 3, [&](uint neighborIdx) {
				const CellData& neighbor = cells[neighborIdx];
				if (!neighbor.hasHit || neighbor.hasClampedHit) {
					isBoundary = true;
					return false;
				}
				return true;
			});

			if (!isBoundary) {
				return;
			}

			const Point3d& hitPoint = cells[idx].hitPoints[0];
			for (const DraftAnglePoint& other : comparisonPoints) {
				if (other.idx == idx) {
					continue;
				}

				if (isWithinPlaneAngle(
						hitPoint,
						other.point,
						normalizedDirection,
						coneCosThreshold,
						eps)) {
					removeCell[idx] = true;
					break;
				}
			}
		});

		parallelFor(allCells, [&](uint idx) {
			if (removeCell[idx]) {
				cells[idx].hasHit = false;
				cells[idx].hasClampedHit = false;
				cells[idx].distance = maxDistance;
				cells[idx].hitPoints.clear();
				removedAny.store(true, std::memory_order_relaxed);
			}
		});
	} while (removedAny.load(std::memory_order_relaxed));
}

static std::vector<CellData> reducePoints(
	std::vector<CellData> cells,
	const GridChoice& grid,
	const vcl::Point3d& direction,
	double draftAngleDegrees,
	float eps,
	double angleThresholdDegrees,
	double maxDistance,
	bool debug = false,
	const std::string& debugResultsSubdir = "",
	vcl::uint* debugStepIndex = nullptr)
{
	using namespace vcl;

	if (cells.size() != grid.rows * grid.cols) {
		return cells;
	}

	std::vector<uint> allCells(cells.size());
	std::iota(allCells.begin(), allCells.end(), 0);
	std::vector<CellData> candidateCells = cells;
	parallelFor(allCells, [&](uint idx) {
		candidateCells[idx].hasHit =
			cells[idx].hasHit && !cells[idx].hasClampedHit;
	});

	erodeHitMaskOnce(candidateCells, grid);
	erodeHitMaskOnce(candidateCells, grid);
	dilateHitMaskOnce(candidateCells, grid);
	dilateHitMaskOnce(candidateCells, grid);

	if (debug && debugStepIndex != nullptr) {
		saveMoldCheckStepMesh( // Step 3
			candidateCells,
			direction,
			debugResultsSubdir,
			*debugStepIndex);
	}
	std::vector<std::vector<uint>> connectedNeighbors =
		removeAngleJumpPoints(
			candidateCells,
			grid,
			direction,
			angleThresholdDegrees,
			eps);

	if (debug && debugStepIndex != nullptr) {
		saveMoldCheckStepMesh( // Step 4
			candidateCells,
			direction,
			debugResultsSubdir,
			*debugStepIndex);
	}

	candidateCells = keepLargestHitComponent(candidateCells, connectedNeighbors);

	if (debug && debugStepIndex != nullptr) {
		saveMoldCheckStepMesh( // Step 5
			candidateCells,
			direction,
			debugResultsSubdir,
			*debugStepIndex);
	}

	removeDraftAngleBoundaryPoints(candidateCells, grid, direction, draftAngleDegrees,eps, maxDistance);

	if (debug && debugStepIndex != nullptr) {
		saveMoldCheckStepMesh( // Step 6
			candidateCells,
			direction,
			debugResultsSubdir,
			*debugStepIndex);
	}

	erodeHitMaskOnce(candidateCells, grid);
	erodeHitMaskOnce(candidateCells, grid);
	erodeHitMaskOnce(candidateCells, grid);
	dilateHitMaskOnce(candidateCells, grid);
	dilateHitMaskOnce(candidateCells, grid);
	dilateHitMaskOnce(candidateCells, grid);

	if (debug && debugStepIndex != nullptr) {
		saveMoldCheckStepMesh( // Step 7
			candidateCells,
			direction,
			debugResultsSubdir,
			*debugStepIndex);
	}

	candidateCells = keepLargestHitComponent(candidateCells, connectedNeighbors);

	if (debug && debugStepIndex != nullptr) {
		saveMoldCheckStepMesh( // Step 8
			candidateCells,
			direction,
			debugResultsSubdir,
			*debugStepIndex);
	}

	candidateCells = cutProtrusions(candidateCells, grid, maxDistance, connectedNeighbors);


	std::vector<CellData> reducedCells = cells;
	parallelFor(allCells, [&](uint idx) {
		if (!candidateCells[idx].hasHit) {
			reducedCells[idx].hasHit = false;
			reducedCells[idx].hasClampedHit = false;
		}
	});

	if (debug && debugStepIndex != nullptr) {
		saveMoldCheckStepMesh( // Step 9
			reducedCells,
			direction,
			debugResultsSubdir,
			*debugStepIndex);
	}

	return reducedCells;
}


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
		if (originalDepthCells[i].hitPoints.empty()) {
			return;
		}

		const Point3d original = originalDepthCells[i].hitPoints[0];

		double requiredT = 0.0;
		bool hasViolation = false;

		for (uint j = 0; j < originalDepthCells.size(); ++j) {
			if (i == j || originalDepthCells[j].hitPoints.empty()) {
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
		depthCells[i].clampedDistance = depthCells[i].distance;
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

	depthCells =
		biharmonicFillWhiteCellsFromRedCells(
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
		3,
		maxDistance);

	if (debugStepIndex != nullptr) {
		saveMoldCheckStepMesh( // Step 11
			depthCells,
			direction,
			debugResultsSubdir,
			*debugStepIndex);
	}

	depthCells =
		biharmonicFillWhiteCellsFromRedCells(
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

	// Only white cells on their upper bound (cyan) remain fixed anchors.
	std::vector<unsigned char> fixedBlueCells(depthCells.size(), false);
	for (uint idx = 0; idx < surfaceCells.size(); ++idx) {
		fixedBlueCells[idx] =
			!surfaceCells[idx].hasHit && depthCells[idx].isMovedForward;
	}

	depthCells = biharmonicFillHitCells(
		surfaceCells,
		depthCells,
		grid,
		direction,
		eps,
		3,
		maxDistance,
		&fixedBlueCells);
	
	if (debugStepIndex != nullptr) {
		saveMoldCheckStepMesh( // Step 13
			depthCells,
			direction,
			debugResultsSubdir,
			*debugStepIndex);
	}
	if (debugStepIndex != nullptr) {
		saveMoldCheckStepMesh( // Step 14
			depthCells,
			direction,
			debugResultsSubdir,
			*debugStepIndex);
	}

	depthCells = fixDepthCellConeViolations(depthCells, direction, coneCosThreshold, eps);

	return depthCells;
}

#endif
