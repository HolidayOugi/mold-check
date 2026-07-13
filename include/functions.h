#ifndef VCL_TEST_EXTERNAL_888_MOLD_CHECK_FUNCTIONS_H
#define VCL_TEST_EXTERNAL_888_MOLD_CHECK_FUNCTIONS_H

#include "helper.h"
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

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>


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

static std::vector<CellData> keepClampedCellsConnectedToCandidates(
	std::vector<CellData> cells,
	const std::vector<CellData>& candidateCells,
	const GridChoice& grid,
	double distanceThreshold)
{
	using namespace vcl;

	std::vector<uint> allCells(cells.size());
	std::iota(allCells.begin(), allCells.end(), 0);
	std::vector<char> keepCells(cells.size(), false);
	std::vector<uint> stack;

	parallelFor(allCells, [&](uint idx) {
		if (candidateCells[idx].hasHit) {
			keepCells[idx] = true;
		}
	});

	for (uint idx = 0; idx < cells.size(); ++idx) {
		if (keepCells[idx]) {
			stack.push_back(idx);
		}
	}

	while (!stack.empty()) {
		const uint idx = stack.back();
		stack.pop_back();

		forEachCrossNeighbor(idx, grid, [&](uint neighborIdx) {
			if (keepCells[neighborIdx]) {
				return true;
			}

			const bool canKeep =
				candidateCells[neighborIdx].hasHit ||
				(cells[neighborIdx].hasHit &&
				 cells[neighborIdx].hasClampedHit);

			if (!canKeep) {
				return true;
			}

			keepCells[neighborIdx] = true;
			stack.push_back(neighborIdx);
			return true;
		});
	}

	std::vector<CellData> filteredCells = cells;
	parallelFor(allCells, [&](uint idx) {
		filteredCells[idx].hasHit = keepCells[idx];
	});

	const std::vector<std::vector<uint>> filteredConnectedNeighbors =
		removeDistanceJumpPoints(
			filteredCells,
			grid,
			distanceThreshold);

	filteredCells =
		keepLargestHitComponent(filteredCells, filteredConnectedNeighbors);

	parallelFor(allCells, [&](uint idx) {
		if (!filteredCells[idx].hasHit) {
			cells[idx].hasHit = false;
			cells[idx].hasClampedHit = false;
		}
	});

	return cells;
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
	double distanceThreshold,
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
		saveMoldCheckStepMesh(
			candidateCells,
			direction,
			debugResultsSubdir,
			*debugStepIndex);
	}

	std::vector<std::vector<uint>> connectedNeighbors =
		removeDistanceJumpPoints(candidateCells,grid, distanceThreshold);

	if (debug && debugStepIndex != nullptr) {
		saveMoldCheckStepMesh(
			candidateCells,
			direction,
			debugResultsSubdir,
			*debugStepIndex);
	}

	candidateCells = keepLargestHitComponent(candidateCells, connectedNeighbors);

	if (debug && debugStepIndex != nullptr) {
		saveMoldCheckStepMesh(
			candidateCells,
			direction,
			debugResultsSubdir,
			*debugStepIndex);
	}

	removeDraftAngleBoundaryPoints(candidateCells, grid, direction, draftAngleDegrees,eps, maxDistance);

	if (debug && debugStepIndex != nullptr) {
		saveMoldCheckStepMesh(
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
		saveMoldCheckStepMesh(
			candidateCells,
			direction,
			debugResultsSubdir,
			*debugStepIndex);
	}

	candidateCells = keepLargestHitComponent(candidateCells, connectedNeighbors);

	if (debug && debugStepIndex != nullptr) {
		saveMoldCheckStepMesh(
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
		saveMoldCheckStepMesh(
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

static std::vector<CellData> biharmonicFillHitCells(
	std::vector<CellData> cells,
	std::vector<CellData> depthCells,
	const GridChoice& grid,
	const vcl::Point3d& direction,
	double eps,
	vcl::uint collarRadius = 3,
	double maxDistance = std::numeric_limits<double>::infinity())
{
	using namespace vcl;

	if (depthCells.size() != grid.rows * grid.cols) {
		return depthCells;
	}

	const auto nearestHitCellDistance = [&](uint idx, uint radius) {
		const uint row = idx / grid.cols;
		const uint col = idx % grid.cols;
		uint nearestDistance = radius + 1;

		const int minRow =
			std::max<int>(0, static_cast<int>(row) - static_cast<int>(radius));
		const int maxRow =
			std::min<int>(
				static_cast<int>(grid.rows) - 1,
				static_cast<int>(row) + static_cast<int>(radius));
		const int minCol =
			std::max<int>(0, static_cast<int>(col) - static_cast<int>(radius));
		const int maxCol =
			std::min<int>(
				static_cast<int>(grid.cols) - 1,
				static_cast<int>(col) + static_cast<int>(radius));

		for (int r = minRow; r <= maxRow; ++r) {
			for (int c = minCol; c <= maxCol; ++c) {
				const uint neighborIdx =
					static_cast<uint>(r) * grid.cols +
					static_cast<uint>(c);
				if (depthCells[neighborIdx].hasHit) {
					const uint rowDistance = static_cast<uint>(
						std::abs(r - static_cast<int>(row)));
					const uint colDistance = static_cast<uint>(
						std::abs(c - static_cast<int>(col)));
					nearestDistance = std::min(
						nearestDistance,
						std::max(rowDistance, colDistance));
				}
			}
		}

		return nearestDistance;
	};

	std::vector<int> unknownVarIds(depthCells.size(), -1);
	std::vector<char> isFixedCollar(depthCells.size(), false);
	std::vector<double> originalDistanceWeights(depthCells.size(), 1.0);
	std::vector<uint> unknownIds;
	std::vector<uint> fixedIds;
	unknownIds.reserve(depthCells.size());

	for (uint idx = 0; idx < depthCells.size(); ++idx) {
		if (!std::isfinite(depthCells[idx].distance)) {
			continue;
		}

		const uint hitDistance =
			depthCells[idx].hasHit ?
				0 :
				nearestHitCellDistance(idx, collarRadius);

		const bool isWeightedCollar =
			!depthCells[idx].hasHit &&
			collarRadius > 0 &&
			hitDistance <= collarRadius;

		if (depthCells[idx].hasHit) {
			originalDistanceWeights[idx] = 0.0;
		}
		else if (isWeightedCollar) {
			originalDistanceWeights[idx] =
				static_cast<double>(hitDistance) /
				static_cast<double>(collarRadius + 1);
		}

		if (depthCells[idx].hasClampedHit ||
			(!depthCells[idx].hasHit && !isWeightedCollar)) {
			continue;
		}

		unknownVarIds[idx] = static_cast<int>(unknownIds.size());
		unknownIds.push_back(idx);
	}

	if (unknownIds.empty()) {
		return depthCells;
	}

	for (uint idx = 0; idx < depthCells.size(); ++idx) {
		if (unknownVarIds[idx] >= 0 ||
			!std::isfinite(depthCells[idx].distance)) {
			continue;
		}

		const uint row = idx / grid.cols;
		const uint col = idx % grid.cols;
		bool touchesUnknown = false;

		if (collarRadius == 0) {
			const auto checkNeighbor = [&](uint neighborIdx) {
				if (unknownVarIds[neighborIdx] >= 0) {
					touchesUnknown = true;
				}
			};

			if (col > 0) {
				checkNeighbor(idx - 1);
			}
			if (!touchesUnknown && col + 1 < grid.cols) {
				checkNeighbor(idx + 1);
			}
			if (!touchesUnknown && row > 0) {
				checkNeighbor(idx - grid.cols);
			}
			if (!touchesUnknown && row + 1 < grid.rows) {
				checkNeighbor(idx + grid.cols);
			}
		}
		else {
			const int minRow = std::max<int>(
				0,
				static_cast<int>(row) - static_cast<int>(collarRadius));
			const int maxRow = std::min<int>(
				static_cast<int>(grid.rows) - 1,
				static_cast<int>(row) + static_cast<int>(collarRadius));
			const int minCol = std::max<int>(
				0,
				static_cast<int>(col) - static_cast<int>(collarRadius));
			const int maxCol = std::min<int>(
				static_cast<int>(grid.cols) - 1,
				static_cast<int>(col) + static_cast<int>(collarRadius));

			for (int r = minRow; r <= maxRow && !touchesUnknown; ++r) {
				for (int c = minCol; c <= maxCol; ++c) {
					const uint neighborIdx =
						static_cast<uint>(r) * grid.cols +
						static_cast<uint>(c);
					if (unknownVarIds[neighborIdx] >= 0) {
						touchesUnknown = true;
						break;
					}
				}
			}
		}

		if (touchesUnknown) {
			isFixedCollar[idx] = true;
			fixedIds.push_back(idx);
		}
	}

	if (fixedIds.empty()) {
		return depthCells;
	}

	std::cout << "  biharmonic unknown cells: " << unknownIds.size()
			  << ", fixed anchor cells: " << fixedIds.size() << "\n";
	std::cout.flush();

	std::vector<uint> laplacianRowCellIds = unknownIds;
	laplacianRowCellIds.insert(
		laplacianRowCellIds.end(),
		fixedIds.begin(),
		fixedIds.end());

	const Eigen::Index rowCount =
		static_cast<Eigen::Index>(laplacianRowCellIds.size());
	const Eigen::Index unknownCount =
		static_cast<Eigen::Index>(unknownIds.size());

	std::vector<Eigen::Triplet<double>> laplacianTriplets;
	laplacianTriplets.reserve(laplacianRowCellIds.size() * 5);
	Eigen::VectorXd fixedRhs = Eigen::VectorXd::Zero(rowCount);

	for (uint rowIdx = 0; rowIdx < laplacianRowCellIds.size(); ++rowIdx) {
		const uint cellIdx = laplacianRowCellIds[rowIdx];
		const uint cellRow = cellIdx / grid.cols;
		const uint cellCol = cellIdx % grid.cols;
		uint usedNeighborCount = 0;

		const auto addNeighbor = [&](uint neighborIdx) {
			if (unknownVarIds[neighborIdx] >= 0) {
				laplacianTriplets.emplace_back(
					static_cast<int>(rowIdx),
					unknownVarIds[neighborIdx],
					1.0);
				++usedNeighborCount;
				return;
			}

			if (isFixedCollar[neighborIdx]) {
				fixedRhs(static_cast<Eigen::Index>(rowIdx)) -=
					depthCells[neighborIdx].distance;
				++usedNeighborCount;
			}
		};

		if (cellCol > 0) {
			addNeighbor(cellIdx - 1);
		}
		if (cellCol + 1 < grid.cols) {
			addNeighbor(cellIdx + 1);
		}
		if (cellRow > 0) {
			addNeighbor(cellIdx - grid.cols);
		}
		if (cellRow + 1 < grid.rows) {
			addNeighbor(cellIdx + grid.cols);
		}

		if (usedNeighborCount == 0) {
			if (unknownVarIds[cellIdx] >= 0) {
				laplacianTriplets.emplace_back(
					static_cast<int>(rowIdx),
					unknownVarIds[cellIdx],
					1.0);
				fixedRhs(static_cast<Eigen::Index>(rowIdx)) =
					depthCells[cellIdx].distance;
			}
			continue;
		}

		if (unknownVarIds[cellIdx] >= 0) {
			laplacianTriplets.emplace_back(
				static_cast<int>(rowIdx),
				unknownVarIds[cellIdx],
				-static_cast<double>(usedNeighborCount));
		}
		else {
			fixedRhs(static_cast<Eigen::Index>(rowIdx)) +=
				static_cast<double>(usedNeighborCount) *
				depthCells[cellIdx].distance;
		}
	}

	Eigen::SparseMatrix<double> laplacian(rowCount, unknownCount);
	laplacian.setFromTriplets(
		laplacianTriplets.begin(),
		laplacianTriplets.end());

	Eigen::SparseMatrix<double> system =
		laplacian.transpose() * laplacian;
	Eigen::VectorXd rhs =
		laplacian.transpose() * fixedRhs;

	const double regularization =
		std::max(1e-12, eps * eps);
	for (uint i = 0; i < unknownIds.size(); ++i) {
		system.coeffRef(
			static_cast<Eigen::Index>(i),
			static_cast<Eigen::Index>(i)) += regularization;
		rhs(static_cast<Eigen::Index>(i)) +=
			regularization * depthCells[unknownIds[i]].distance;
	}
	system.makeCompressed();

	std::cout << "  biharmonic sparse solve start\n";
	std::cout.flush();

	Eigen::ConjugateGradient<
		Eigen::SparseMatrix<double>,
		Eigen::Lower | Eigen::Upper> solver;
	solver.setTolerance(1e-8);
	solver.setMaxIterations(
		static_cast<int>(std::max<size_t>(1000, unknownIds.size() * 2)));
	solver.compute(system);

	if (solver.info() != Eigen::Success) {
		return depthCells;
	}

	const Eigen::VectorXd solvedDistances = solver.solve(rhs);

	std::cout << "  biharmonic sparse solve done. Iterations: "
			  << solver.iterations()
			  << ", error: " << solver.error() << "\n";
	std::cout.flush();

	if (solver.info() != Eigen::Success) {
		return depthCells;
	}

	for (uint i = 0; i < unknownIds.size(); ++i) {
		double distance =
			solvedDistances(static_cast<Eigen::Index>(i));

		if (!std::isfinite(distance)) {
			continue;
		}

		bool ignoreMax = false;

		CellData& cell = depthCells[unknownIds[i]];
		const double originalDistance = cell.distance;

		const std::vector<uint> neighbors = squareNeighborIndices(unknownIds[i], grid, 5);
		for (uint neighborIdx : neighbors) {
			if (!depthCells[neighborIdx].hasHit) {
				ignoreMax = true;
				break;
			}
		}

		//ignoreMax = true; //debug REMOVE!!!!

		if (!ignoreMax) {
			cell.distance = std::max(
				distance,
				cells[unknownIds[i]].distance + 0.003 * maxDistance);
		}
		else {
			cell.distance = distance;
		}

		const double originalDistanceWeight =
			originalDistanceWeights[unknownIds[i]];
		if (originalDistanceWeight > 0.0) {
			cell.distance =
				originalDistanceWeight * originalDistance +
				(1.0 - originalDistanceWeight) * cell.distance;
		}

		cell.hitPoints = {cell.cellCenter + direction * cell.distance};

		if (cell.hasClampedHit) {
			cell.clampedDistance = cell.distance;
		}
	}

	std::cout << "  biharmonic done\n";
	std::cout.flush();

	return depthCells;
}

static std::vector<double> pushPull(
	const std::vector<CellData>& cells,
	const GridChoice& grid)
{
	using namespace vcl;

	std::vector<PullPushLevel> pyramid;

	PullPushLevel base;
	base.rows = grid.rows;
	base.cols = grid.cols;
	base.distances.resize(cells.size(), 0.0);
	base.weights.resize(cells.size(), 0.0);

	std::vector<uint> baseIndices(cells.size());
	std::iota(baseIndices.begin(), baseIndices.end(), 0);

	parallelFor(baseIndices, [&](uint i) {
		if (cells[i].hasHit && !cells[i].hasClampedHit) {
			base.distances[i] = cells[i].clampedDistance;
			base.weights[i] = 1.0;
		}
	});

	const uint fixedHitCount = static_cast<uint>(
		std::count_if(cells.begin(), cells.end(), [](const CellData& cell) {
			return cell.hasHit && !cell.hasClampedHit;
		}));

	if (fixedHitCount == 0) {
		return base.distances;
	}

	pyramid.push_back(std::move(base));

	while (pyramid.back().rows > 1 || pyramid.back().cols > 1) {
		const PullPushLevel& fine = pyramid.back();

		PullPushLevel coarse;
		coarse.rows = (fine.rows + 1) / 2;
		coarse.cols = (fine.cols + 1) / 2;
		coarse.distances.resize(coarse.rows * coarse.cols, 0.0);
		coarse.weights.resize(coarse.rows * coarse.cols, 0.0);

		std::vector<uint> coarseIndices(coarse.rows * coarse.cols);
		std::iota(coarseIndices.begin(), coarseIndices.end(), 0);

		parallelFor(coarseIndices, [&](uint coarseIndex) {
			const uint row = coarseIndex / coarse.cols;
			const uint col = coarseIndex % coarse.cols;
			double weightedDistanceSum = 0.0;
			double weightSum = 0.0;

			for (uint rowOffset = 0; rowOffset < 2; ++rowOffset) {
				for (uint colOffset = 0; colOffset < 2; ++colOffset) {
					const uint fineRow = row * 2 + rowOffset;
					const uint fineCol = col * 2 + colOffset;

					if (fineRow >= fine.rows || fineCol >= fine.cols) {
						continue;
					}

					const uint fineIndex = fineRow * fine.cols + fineCol;
					const double weight = fine.weights[fineIndex];

					if (weight <= 0.0) {
						continue;
					}

					weightedDistanceSum += fine.distances[fineIndex] * weight;
					weightSum += weight;
				}
			}

			if (weightSum > 0.0) {
				coarse.distances[coarseIndex] =
					weightedDistanceSum / weightSum;
				coarse.weights[coarseIndex] = weightSum;
			}
		});

		pyramid.push_back(std::move(coarse));
	}

	for (uint level = static_cast<uint>(pyramid.size() - 1);
		 level > 0;
		 --level) {
		const PullPushLevel& coarse = pyramid[level];
		PullPushLevel& fine = pyramid[level - 1];

		std::vector<uint> fineIndices(fine.rows * fine.cols);
		std::iota(fineIndices.begin(), fineIndices.end(), 0);

		parallelFor(fineIndices, [&](uint index) {
			const uint row = index / fine.cols;
			const uint col = index % fine.cols;

			if (level == 1 && cells[index].hasHit && !cells[index].hasClampedHit) {
				fine.distances[index] = cells[index].clampedDistance;
				fine.weights[index] = 1.0;
				return;
			}

			if (level > 1 && fine.weights[index] > 0.0) {
				return;
			}

			fine.distances[index] =
				interpolateFromCoarseLevel(
					coarse,
					row,
					col,
					fine.rows,
					fine.cols);

			if (level == 1 && cells[index].hasClampedHit) {
				fine.distances[index] =
					std::min(fine.distances[index], cells[index].clampedDistance);
			}

			fine.weights[index] = 1.0;
		});
	}

	return pyramid[0].distances;
}

static std::vector<CellData> successiveOverRelaxation(
	std::vector<CellData> depthCells,
	const std::vector<CellData>& cells,
	const vcl::Point3d& direction,
	const GridChoice& grid,
	vcl::uint maxIterations,
	double omega,
	double eps)
{
	using namespace vcl;

	if (depthCells.size() != cells.size() ||
		depthCells.size() != grid.rows * grid.cols ||
		omega <= 0.0 ||
		omega >= 2.0) {
		return depthCells;
	}

	std::vector<uint> allCells(depthCells.size());
	std::iota(allCells.begin(), allCells.end(), 0);
	std::vector<double> changes(depthCells.size(), 0.0);

	for (uint iteration = 0; iteration < maxIterations; ++iteration) {
		std::fill(changes.begin(), changes.end(), 0.0);

		for (uint color = 0; color < 2; ++color) {
			parallelFor(allCells, [&](uint idx) {
				if (cells[idx].hasHit && !cells[idx].hasClampedHit) {
					return;
				}

				const uint row = idx / grid.cols;
				const uint col = idx % grid.cols;
				if ((row + col) % 2 != color) {
					return;
				}

				const std::vector<uint> neighbors = squareNeighborIndices(idx, grid, 3);
				if (neighbors.empty()) {
					return;
				}

				double neighborDistanceSum = 0.0;
				for (uint neighborIdx : neighbors) {
					neighborDistanceSum += depthCells[neighborIdx].distance;
				}

				const double gaussSeidelDistance =
					neighborDistanceSum / static_cast<double>(neighbors.size());
				const double oldDistance = depthCells[idx].distance;
				const double relaxedDistance =
					(1.0 - omega) * oldDistance +
					omega * gaussSeidelDistance;

				depthCells[idx].distance = relaxedDistance;
				changes[idx] = std::abs(relaxedDistance - oldDistance);
			});
		}

		const double maxChange =
			*std::max_element(changes.begin(), changes.end());
		if (maxChange < eps) {
			break;
		}
	}

	parallelFor(allCells, [&](uint i) {
		depthCells[i].hitPoints = {
			depthCells[i].cellCenter + direction * depthCells[i].distance};
	});

	return depthCells;
}

static std::vector<CellData> smoothMeshEdges(
	std::vector<CellData> depthCells,
	const vcl::Point3d& direction,
	const GridChoice& grid)
{
	using namespace vcl;

	if (depthCells.size() != grid.rows * grid.cols) {
		return depthCells;
	}

	const std::vector<CellData> originalDepthCells = depthCells;
	std::vector<uint> allCells(depthCells.size());
	std::iota(allCells.begin(), allCells.end(), 0);

	parallelFor(allCells, [&](uint idx) {
		if (!originalDepthCells[idx].hasHit) {
			return;
		}

		double clampedDistance = originalDepthCells[idx].distance;
		bool shouldClamp = false;

		for (uint neighborIdx : squareNeighborIndices(idx, grid, 3)) {
			if (originalDepthCells[neighborIdx].hasHit) {
				continue;
			}

			if (clampedDistance > originalDepthCells[neighborIdx].distance) {
				clampedDistance = originalDepthCells[neighborIdx].distance;
				shouldClamp = true;
			}
		}

		if (!shouldClamp) {
			return;
		}

		depthCells[idx].distance = clampedDistance;
		depthCells[idx].clampedDistance = clampedDistance;
		depthCells[idx].hitPoints = {
			depthCells[idx].cellCenter + direction * clampedDistance};
		depthCells[idx].hasClampedHit = true;
	});

	return depthCells;
}

static std::vector<CellData> smoothMeshPits(
	std::vector<CellData> depthCells,
	const vcl::Point3d& direction,
	const GridChoice& grid)
{
	using namespace vcl;

	if (depthCells.size() != grid.rows * grid.cols) {
		return depthCells;
	}

	const std::vector<CellData> originalDepthCells = depthCells;
	std::vector<uint> allCells(depthCells.size());
	std::iota(allCells.begin(), allCells.end(), 0);

	parallelFor(allCells, [&](uint idx) {
		if (!originalDepthCells[idx].hasHit ||
			originalDepthCells[idx].hasClampedHit) {
			return;
		}

		const std::vector<uint> neighbors =
			squareNeighborIndices(idx, grid, 3);

		if (neighbors.empty()) {
			return;
		}

		uint pitNeighborCount = 0;
		double maxClampedDistance =
			-std::numeric_limits<double>::infinity();

		for (uint neighborIdx : neighbors) {
			if (!originalDepthCells[neighborIdx].hasHit ||
				originalDepthCells[neighborIdx].hasClampedHit) {
				++pitNeighborCount;
			}

			if (!originalDepthCells[neighborIdx].hasClampedHit) {
				continue;
			}

			maxClampedDistance =
				std::max(
					maxClampedDistance,
					originalDepthCells[neighborIdx].distance);
		}

		const bool hasClampedNeighbor =
			maxClampedDistance > -std::numeric_limits<double>::infinity();

		const bool mostlySurroundedByPitNeighbors =
			pitNeighborCount + 1 > neighbors.size() / 2;

		if (!hasClampedNeighbor ||
			!mostlySurroundedByPitNeighbors ||
			originalDepthCells[idx].distance <= maxClampedDistance) {
			return;
		}

		depthCells[idx].distance = maxClampedDistance;
		depthCells[idx].clampedDistance = maxClampedDistance;
		depthCells[idx].hitPoints = {
			depthCells[idx].cellCenter + direction * maxClampedDistance};
		depthCells[idx].hasClampedHit = true;
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

	std::vector<CellData> depthCells = cells;

	if (std::none_of(cells.begin(), cells.end(), [](const CellData& cell) {
			return cell.hasHit;
		})) {
		return depthCells;
	}

	const std::vector<double> distances = pushPull(cells, grid);

	for (uint i = 0; i < depthCells.size(); ++i) {
		depthCells[i].distance = distances[i];

		if (cells[i].hasClampedHit) {
			depthCells[i].distance =
				std::min(depthCells[i].distance, cells[i].clampedDistance);
		}

		depthCells[i].hitPoints = {
			depthCells[i].cellCenter + direction * depthCells[i].distance
		};
		depthCells[i].hasHit = cells[i].hasHit;
	}

	if (debugStepIndex != nullptr) {
		saveMoldCheckStepMesh(
			depthCells,
			direction,
			debugResultsSubdir,
			*debugStepIndex);
	}

	depthCells = successiveOverRelaxation(depthCells, cells, direction, grid, 1000, 1.6, eps);

	if (debugStepIndex != nullptr) {
		saveMoldCheckStepMesh(
			depthCells,
			direction,
			debugResultsSubdir,
			*debugStepIndex);
	}

	depthCells = fixDepthCellConeViolations(depthCells, direction, coneCosThreshold, eps);

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

	const std::filesystem::path debugOutputDir =
			std::filesystem::path(RESULTS_PATH) /
			debugResultsSubdir;

	std::filesystem::create_directories(debugOutputDir);
	
	const std::string base =
			(debugOutputDir / "mold_check").string();
	saveMesh(depthPointsMesh, base + "_original_mold_points.ply");

	depthCells = biharmonicFillHitCells(
		cells,
		depthCells,
		grid,
		direction,
		eps,
		3,
		maxDistance);
	
	if (debugStepIndex != nullptr) {
		saveMoldCheckStepMesh(
			depthCells,
			direction,
			debugResultsSubdir,
			*debugStepIndex);
	}

	depthCells = fixDepthCellConeViolations(depthCells, direction, coneCosThreshold, eps);

	return depthCells;
}

#endif
