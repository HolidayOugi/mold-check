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
#include <queue>
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

static std::vector<CellData> biharmonicFillWhiteCellsFromRedCells(
	const std::vector<CellData>& cells,
	std::vector<CellData> depthCells,
	const GridChoice& grid,
	const vcl::Point3d& direction,
	double eps,
	double maxDistance = std::numeric_limits<double>::infinity())
{
	using namespace vcl;

	if (depthCells.size() != cells.size() ||
		depthCells.size() != grid.rows * grid.cols) {
		return depthCells;
	}

	std::vector<int> unknownVarIds(depthCells.size(), -1);
	std::vector<char> isFixedRed(depthCells.size(), false);
	std::vector<uint> unknownIds;
	std::vector<uint> fixedIds;
	unknownIds.reserve(depthCells.size());

	for (uint idx = 0; idx < depthCells.size(); ++idx) {
		if (!std::isfinite(depthCells[idx].distance)) {
			continue;
		}

		if (cells[idx].hasHit && !cells[idx].hasClampedHit) {
			isFixedRed[idx] = true;
			continue;
		}

		if (!cells[idx].hasHit) {
			unknownVarIds[idx] = static_cast<int>(unknownIds.size());
			unknownIds.push_back(idx);
		}
	}

	if (unknownIds.empty()) {
		return depthCells;
	}

	for (uint idx = 0; idx < depthCells.size(); ++idx) {
		if (!isFixedRed[idx]) {
			continue;
		}

		bool touchesUnknown = false;
		forEachCrossNeighbor(idx, grid, [&](uint neighborIdx) {
			if (unknownVarIds[neighborIdx] >= 0) {
				touchesUnknown = true;
				return false;
			}
			return true;
		});

		if (touchesUnknown) {
			fixedIds.push_back(idx);
		}
	}

	if (fixedIds.empty()) {
		return depthCells;
	}

	std::cout << "  biharmonic white cells from red cells. Unknown cells: "
			  << unknownIds.size()
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

			if (isFixedRed[neighborIdx]) {
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

	Eigen::VectorXd solvedDistances(unknownCount);
	int solveIterations = 0;
	double solveError = std::numeric_limits<double>::infinity();

	if (!std::isfinite(maxDistance)) {
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

		solvedDistances = solver.solve(rhs);
		solveIterations = solver.iterations();
		solveError = solver.error();

		if (solver.info() != Eigen::Success) {
			return depthCells;
		}
	}
	else {
		Eigen::VectorXd upperBounds(unknownCount);
		for (uint i = 0; i < unknownIds.size(); ++i) {
			const uint cellIdx = unknownIds[i];
			const double surfaceDistance = cells[cellIdx].distance;
			const double upperBound =
				std::isfinite(surfaceDistance) ?
					surfaceDistance - 0.01 * maxDistance :
					std::numeric_limits<double>::infinity();
			upperBounds(static_cast<Eigen::Index>(i)) = upperBound;
			solvedDistances(static_cast<Eigen::Index>(i)) = std::min(
				depthCells[cellIdx].distance,
				upperBound);
		}

		double lipschitzConstant = 0.0;
		for (Eigen::Index outer = 0; outer < system.outerSize(); ++outer) {
			double absoluteColumnSum = 0.0;
			for (Eigen::SparseMatrix<double>::InnerIterator entry(system, outer);
				 entry;
				 ++entry) {
				absoluteColumnSum += std::abs(entry.value());
			}
			lipschitzConstant = std::max(
				lipschitzConstant,
				absoluteColumnSum);
		}

		if (!std::isfinite(lipschitzConstant) ||
			lipschitzConstant <= 0.0 ||
			!solvedDistances.allFinite()) {
			return depthCells;
		}

		Eigen::VectorXd acceleratedDistances = solvedDistances;
		double momentum = 1.0;
		const double inverseLipschitz = 1.0 / lipschitzConstant;
		const size_t requestedIterations =
			std::max<size_t>(1000, unknownIds.size() * 2);
		const int maximumIterations = static_cast<int>(
			std::min<size_t>(20000, requestedIterations));

		for (; solveIterations < maximumIterations; ++solveIterations) {
			const Eigen::VectorXd gradient =
				system * acceleratedDistances - rhs;
			Eigen::VectorXd nextDistances =
				(acceleratedDistances - inverseLipschitz * gradient)
					.cwiseMin(upperBounds);

			if (!nextDistances.allFinite()) {
				return depthCells;
			}

			solveError =
				(nextDistances - solvedDistances).norm() /
				std::max(1.0, solvedDistances.norm());
			if (solveError <= 1e-8) {
				solvedDistances = nextDistances;
				++solveIterations;
				break;
			}

			const double nextMomentum =
				0.5 * (1.0 + std::sqrt(1.0 + 4.0 * momentum * momentum));
			acceleratedDistances =
				nextDistances +
				((momentum - 1.0) / nextMomentum) *
					(nextDistances - solvedDistances);
			solvedDistances = nextDistances;
			momentum = nextMomentum;
		}

		const Eigen::VectorXd projectedDistances =
			(solvedDistances - inverseLipschitz *
				(system * solvedDistances - rhs))
				.cwiseMin(upperBounds);
		solveError =
			(projectedDistances - solvedDistances).norm() /
			std::max(1.0, solvedDistances.norm());
	}

	std::cout << "  biharmonic sparse solve done. Iterations: "
			  << solveIterations
			  << ", error: " << solveError << "\n";
	std::cout.flush();

	for (uint i = 0; i < unknownIds.size(); ++i) {
		const uint cellIdx = unknownIds[i];
		const double distance =
			solvedDistances(static_cast<Eigen::Index>(i));

		if (!std::isfinite(distance)) {
			continue;
		}

		CellData& cell = depthCells[cellIdx];
		cell.distance = distance;
		cell.hitPoints = {cell.cellCenter + direction * cell.distance};

		if (std::isfinite(maxDistance) &&
			std::isfinite(cells[cellIdx].distance)) {
			const double upperBound =
				cells[cellIdx].distance - 0.01 * maxDistance;
			const double activeTolerance = std::max(
				1e-10,
				10.0 * static_cast<double>(eps) *
					std::max(1.0, std::abs(upperBound)));
			cell.isMovedForward =
				distance >= upperBound - activeTolerance;
		}
	}

	std::cout << "  biharmonic done\n";
	std::cout.flush();

	return depthCells;
}

// Minimizes the squared grid Laplacian over orange hit cells and the nearby
// white collar. The first call keeps the original unconstrained solve; the
// final call jointly applies outside bounds to whites and inside bounds to
// orange cells from the third layer onward.
static std::vector<CellData> biharmonicFillHitCells(
	std::vector<CellData> cells,
	std::vector<CellData> depthCells,
	const GridChoice& grid,
	const vcl::Point3d& direction,
	double eps,
	vcl::uint collarRadius = 3,
	double maxDistance = std::numeric_limits<double>::infinity(),
	const std::vector<unsigned char>* fixedWhiteAnchors = nullptr)
{
	using namespace vcl;

	if (depthCells.size() != cells.size() ||
		depthCells.size() != grid.rows * grid.cols) {
		return depthCells;
	}

	// Returns the Chebyshev distance from a white cell to the hit region.
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

	// Detects orange cells whose first white neighbor is at least N layers away.
	const auto isAtLeastDistanceFromWhiteBoundary = [&](
		uint idx,
		uint minimumDistance) {
		if (!depthCells[idx].hasHit || minimumDistance == 0) {
			return false;
		}

		const uint row = idx / grid.cols;
		const uint col = idx % grid.cols;
		const int searchRadius = static_cast<int>(minimumDistance) - 1;
		const int minRow = std::max<int>(0, static_cast<int>(row) - searchRadius);
		const int maxRow = std::min<int>(
			static_cast<int>(grid.rows) - 1,
			static_cast<int>(row) + searchRadius);
		const int minCol = std::max<int>(0, static_cast<int>(col) - searchRadius);
		const int maxCol = std::min<int>(
			static_cast<int>(grid.cols) - 1,
			static_cast<int>(col) + searchRadius);

		for (int r = minRow; r <= maxRow; ++r) {
			for (int c = minCol; c <= maxCol; ++c) {
				const uint neighborIdx =
					static_cast<uint>(r) * grid.cols + static_cast<uint>(c);
				if (!depthCells[neighborIdx].hasHit) {
					return false;
				}
			}
		}

		return true;
	};

	// Split cells into solver variables and fixed boundary anchors.
	std::vector<int> unknownVarIds(depthCells.size(), -1);
	std::vector<char> isFixedCollar(depthCells.size(), false);
	std::vector<double> originalDistanceWeights(depthCells.size(), 1.0);
	std::vector<uint> unknownIds;
	std::vector<uint> fixedIds;
	unknownIds.reserve(depthCells.size());

	const double depthConstraintMargin =
		std::isfinite(maxDistance) ?
			0.003 * maxDistance :
			std::max<double>(1e-6, 100.0 * eps);

	// A valid anchor mask identifies the final box-constrained solve.
	const bool useBoxConstraints =
		fixedWhiteAnchors != nullptr &&
		fixedWhiteAnchors->size() == depthCells.size() &&
		std::isfinite(maxDistance);

	// Detect complete depressions in the incoming cutting surface. Distances
	// grow away from the ray plane, so a geometric depression is a local basin
	// of the cutting surface but a local maximum in this scalar field.
	std::vector<double> basinUpperBounds(
		depthCells.size(),
		std::numeric_limits<double>::infinity());

	if (useBoxConstraints && grid.rows > 1 && grid.cols > 1) {
		std::vector<double> filledDistances(depthCells.size());
		std::vector<unsigned char> floodVisited(depthCells.size(), false);
		std::priority_queue<std::pair<double, uint>> floodFront;

		for (uint idx = 0; idx < depthCells.size(); ++idx) {
			filledDistances[idx] = depthCells[idx].distance;
		}

		const auto addFloodBoundary = [&](uint idx) {
			if (floodVisited[idx]) {
				return;
			}
			floodVisited[idx] = true;
			if (std::isfinite(filledDistances[idx])) {
				floodFront.push({filledDistances[idx], idx});
			}
		};

		for (uint row = 0; row < grid.rows; ++row) {
			addFloodBoundary(row * grid.cols);
			addFloodBoundary(row * grid.cols + grid.cols - 1);
		}
		for (uint col = 0; col < grid.cols; ++col) {
			addFloodBoundary(col);
			addFloodBoundary((grid.rows - 1) * grid.cols + col);
		}

		const auto visitFloodNeighbor = [&](
			uint neighborIdx,
			double spillDistance) {
			if (floodVisited[neighborIdx]) {
				return;
			}
			floodVisited[neighborIdx] = true;
			if (!std::isfinite(depthCells[neighborIdx].distance)) {
				return;
			}

			filledDistances[neighborIdx] = std::min(
				depthCells[neighborIdx].distance,
				spillDistance);
			floodFront.push({filledDistances[neighborIdx], neighborIdx});
		};

		while (!floodFront.empty()) {
			const auto [spillDistance, idx] = floodFront.top();
			floodFront.pop();

			const uint row = idx / grid.cols;
			const uint col = idx % grid.cols;
			if (col > 0) {
				visitFloodNeighbor(idx - 1, spillDistance);
			}
			if (col + 1 < grid.cols) {
				visitFloodNeighbor(idx + 1, spillDistance);
			}
			if (row > 0) {
				visitFloodNeighbor(idx - grid.cols, spillDistance);
			}
			if (row + 1 < grid.rows) {
				visitFloodNeighbor(idx + grid.cols, spillDistance);
			}
		}

		const double minimumCellSide = std::min(grid.sideU, grid.sideV);
		const double minimumBasinProminence = std::max(
			0.5 * minimumCellSide,
			0.003 * maxDistance);
		const double basinTolerance = std::max(
			1e-10,
			10.0 * eps);

		std::vector<unsigned char> componentVisited(
			depthCells.size(), false);
		std::vector<unsigned char> isSelectedBasin(
			depthCells.size(), false);
		uint selectedBasinCount = 0;
		uint selectedBasinCellCount = 0;
		double maximumBasinLift = 0.0;

		for (uint seedIdx = 0; seedIdx < depthCells.size(); ++seedIdx) {
			const double seedLift =
				depthCells[seedIdx].distance - filledDistances[seedIdx];
			if (componentVisited[seedIdx] ||
				!std::isfinite(seedLift) || seedLift <= basinTolerance) {
				continue;
			}

			std::vector<uint> component;
			std::vector<uint> pendingCells = {seedIdx};
			componentVisited[seedIdx] = true;
			double componentMaximumLift = 0.0;
			double componentMaximumOrangeLift = 0.0;

			while (!pendingCells.empty()) {
				const uint idx = pendingCells.back();
				pendingCells.pop_back();
				component.push_back(idx);
				componentMaximumLift = std::max(
					componentMaximumLift,
					depthCells[idx].distance - filledDistances[idx]);

				if (depthCells[idx].hasHit) {
					componentMaximumOrangeLift = std::max(
						componentMaximumOrangeLift,
						depthCells[idx].distance - filledDistances[idx]);
				}

				const uint row = idx / grid.cols;
				const uint col = idx % grid.cols;
				const auto addComponentNeighbor = [&](uint neighborIdx) {
					const double lift = depthCells[neighborIdx].distance -
						filledDistances[neighborIdx];
					if (!componentVisited[neighborIdx] &&
						std::isfinite(lift) && lift > basinTolerance) {
						componentVisited[neighborIdx] = true;
						pendingCells.push_back(neighborIdx);
					}
				};

				if (col > 0) {
					addComponentNeighbor(idx - 1);
				}
				if (col + 1 < grid.cols) {
					addComponentNeighbor(idx + 1);
				}
				if (row > 0) {
					addComponentNeighbor(idx - grid.cols);
				}
				if (row + 1 < grid.rows) {
					addComponentNeighbor(idx + grid.cols);
				}
			}

			if (componentMaximumLift < minimumBasinProminence ||
				componentMaximumOrangeLift < minimumBasinProminence) {
				continue;
			}

			++selectedBasinCount;
			selectedBasinCellCount += static_cast<uint>(component.size());
			maximumBasinLift = std::max(
				maximumBasinLift,
				componentMaximumLift);
			for (uint idx : component) {
				isSelectedBasin[idx] = true;
			}
		}

		// The support collar prevents the final solve from making the filled
		// depression wider by pulling its immediate surroundings deeper.
		std::vector<unsigned char> isBasinSupport = isSelectedBasin;
		for (uint layer = 0; layer < collarRadius; ++layer) {
			std::vector<unsigned char> expandedSupport = isBasinSupport;
			for (uint idx = 0; idx < depthCells.size(); ++idx) {
				if (!isBasinSupport[idx]) {
					continue;
				}
				const uint row = idx / grid.cols;
				const uint col = idx % grid.cols;
				if (col > 0) expandedSupport[idx - 1] = true;
				if (col + 1 < grid.cols) expandedSupport[idx + 1] = true;
				if (row > 0) expandedSupport[idx - grid.cols] = true;
				if (row + 1 < grid.rows) expandedSupport[idx + grid.cols] = true;
			}
			isBasinSupport.swap(expandedSupport);
		}

		for (uint idx = 0; idx < depthCells.size(); ++idx) {
			if (isSelectedBasin[idx]) {
				basinUpperBounds[idx] = filledDistances[idx];
			}
			else if (isBasinSupport[idx]) {
				basinUpperBounds[idx] = depthCells[idx].distance;
			}
		}

		std::cout << "  depth basins constrained: "
				  << selectedBasinCount
				  << ", basin cells: " << selectedBasinCellCount
				  << ", maximum lift: " << maximumBasinLift << "\n";
		std::cout.flush();
	}

	std::vector<unsigned char> isFixedOutsideOrange(depthCells.size(), false);
	for (uint idx = 0; idx < depthCells.size(); ++idx) {
		if (!std::isfinite(depthCells[idx].distance)) {
			continue;
		}

		const bool isFixedWhiteAnchor =
			useBoxConstraints &&
			(*fixedWhiteAnchors)[idx] &&
			!depthCells[idx].hasHit &&
			!std::isfinite(basinUpperBounds[idx]);

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
		// Keep the legacy snap only in the first solve; the final solve uses bounds.
		isFixedOutsideOrange[idx] =
			!useBoxConstraints &&
			depthCells[idx].hasHit &&
			!depthCells[idx].hasClampedHit &&
			std::isfinite(cells[idx].distance) &&
			depthCells[idx].distance < cells[idx].distance &&
			isAtLeastDistanceFromWhiteBoundary(idx, 3);

		if (isFixedOutsideOrange[idx]) {
			CellData& constrainedCell = depthCells[idx];
			double constrainedDistance =
				cells[idx].distance + depthConstraintMargin;

			if (cells[idx].hitPoints.size() > 1) {
				const double lastHitDistance =
					(cells[idx].hitPoints.back() -
					 cells[idx].cellCenter).dot(direction);
				if (std::isfinite(lastHitDistance)) {
					constrainedDistance = std::min(
						constrainedDistance,
						lastHitDistance - 0.003 * maxDistance);
				}
			}

			constrainedCell.distance = constrainedDistance;
			constrainedCell.hitPoints = {
				constrainedCell.cellCenter + direction * constrainedCell.distance};
			constrainedCell.isBiharmonicFilledHit = true;
		}

		if (isFixedWhiteAnchor ||
			isFixedOutsideOrange[idx] ||
			depthCells[idx].hasClampedHit ||
			(!depthCells[idx].hasHit &&
			 !isWeightedCollar &&
			 !std::isfinite(basinUpperBounds[idx]))) {
			continue;
		}

		unknownVarIds[idx] = static_cast<int>(unknownIds.size());
		unknownIds.push_back(idx);
	}

	if (unknownIds.empty()) {
		return depthCells;
	}

	// Retain finite non-variables near the unknown region as Dirichlet anchors.
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

	// Build Laplacian rows for both variables and their fixed boundary.
	std::vector<uint> laplacianRowCellIds = unknownIds;
	laplacianRowCellIds.insert(
		laplacianRowCellIds.end(),
		fixedIds.begin(),
		fixedIds.end());

	const Eigen::Index rowCount =
		static_cast<Eigen::Index>(laplacianRowCellIds.size());
	const Eigen::Index unknownCount =
		static_cast<Eigen::Index>(unknownIds.size());

	// Assemble the four-neighbor grid Laplacian and move anchors to the RHS.
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

	// Squaring the Laplacian produces the biharmonic quadratic energy.
	Eigen::SparseMatrix<double> laplacian(rowCount, unknownCount);
	laplacian.setFromTriplets(
		laplacianTriplets.begin(),
		laplacianTriplets.end());

	Eigen::SparseMatrix<double> system =
		laplacian.transpose() * laplacian;
	Eigen::VectorXd rhs =
		laplacian.transpose() * fixedRhs;

	// A tiny fidelity term removes the null space and keeps the system stable.
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

	Eigen::VectorXd solvedDistances(unknownCount);
	Eigen::VectorXd lowerBounds;
	Eigen::VectorXd upperBounds;
	int solveIterations = 0;
	double solveError = std::numeric_limits<double>::infinity();

	// The first pass uses the original unconstrained conjugate-gradient solve.
	if (!useBoxConstraints) {
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

		solvedDistances = solver.solve(rhs);
		solveIterations = solver.iterations();
		solveError = solver.error();

		if (solver.info() != Eigen::Success) {
			return depthCells;
		}
	}
	else {
		// The final pass minimizes the same energy under per-cell box bounds.
		lowerBounds = Eigen::VectorXd::Constant(
			unknownCount,
			-std::numeric_limits<double>::infinity());
		upperBounds = Eigen::VectorXd::Constant(
			unknownCount,
			std::numeric_limits<double>::infinity());

		for (uint i = 0; i < unknownIds.size(); ++i) {
			const uint cellIdx = unknownIds[i];
			const Eigen::Index variableIdx =
				static_cast<Eigen::Index>(i);
			const CellData& depthCell = depthCells[cellIdx];
			const CellData& surfaceCell = cells[cellIdx];

			if (!depthCell.hasHit) {
				// White collar cells may move but must remain outside the mesh.
				if (std::isfinite(surfaceCell.distance)) {
					upperBounds(variableIdx) =
						surfaceCell.distance - 0.01 * maxDistance;
				}
			}
			else {
				// Layers 1-2 are free; orange cells from layer 3 must stay inside.
				const bool isBasinConstrained =
					std::isfinite(basinUpperBounds[cellIdx]);
				const bool requiresInteriorDepth =
					!isBasinConstrained &&
					isAtLeastDistanceFromWhiteBoundary(cellIdx, 3);
				if (requiresInteriorDepth &&
					std::isfinite(surfaceCell.distance)) {
					lowerBounds(variableIdx) =
						surfaceCell.distance + depthConstraintMargin;
				}

				if (surfaceCell.hitPoints.size() > 1 &&
					std::isfinite(surfaceCell.distance)) {
					const double lastHitDistance =
						(surfaceCell.hitPoints.back() -
						 surfaceCell.cellCenter).dot(direction);
					// The last-hit upper bound prevents exiting the opposite side.
					if (std::isfinite(lastHitDistance) &&
						lastHitDistance > surfaceCell.distance) {
						const double availableMargin = std::min(
							depthConstraintMargin,
							0.5 * (lastHitDistance - surfaceCell.distance));
						upperBounds(variableIdx) =
							lastHitDistance - availableMargin;
						if (requiresInteriorDepth) {
							lowerBounds(variableIdx) =
								surfaceCell.distance + availableMargin;
						}
					}
				}
			}

			// The basin ceiling is part of this solve. Basin cells rise to at
			// least their spill surface; support cells may rise but cannot deepen.
			if (std::isfinite(basinUpperBounds[cellIdx])) {
				upperBounds(variableIdx) = std::min(
					upperBounds(variableIdx),
					basinUpperBounds[cellIdx]);
				if (lowerBounds(variableIdx) > upperBounds(variableIdx)) {
					lowerBounds(variableIdx) =
					-std::numeric_limits<double>::infinity();
				}
			}

			// Start from a feasible version of the incoming depth.
			double initialDistance = depthCell.distance;
			initialDistance = std::max(
				initialDistance,
				lowerBounds(variableIdx));
			initialDistance = std::min(
				initialDistance,
				upperBounds(variableIdx));
			solvedDistances(variableIdx) = initialDistance;
		}

		// This row-sum bound gives a safe step size for projected FISTA.
		double lipschitzConstant = 0.0;
		for (Eigen::Index outer = 0; outer < system.outerSize(); ++outer) {
			double absoluteColumnSum = 0.0;
			for (Eigen::SparseMatrix<double>::InnerIterator entry(system, outer);
				 entry;
				 ++entry) {
				absoluteColumnSum += std::abs(entry.value());
			}
			lipschitzConstant = std::max(
				lipschitzConstant,
				absoluteColumnSum);
		}

		if (!std::isfinite(lipschitzConstant) ||
			lipschitzConstant <= 0.0 ||
			!solvedDistances.allFinite()) {
			return depthCells;
		}

		Eigen::VectorXd acceleratedDistances = solvedDistances;
		double momentum = 1.0;
		const double inverseLipschitz = 1.0 / lipschitzConstant;
		const size_t requestedIterations =
			std::max<size_t>(1000, unknownIds.size() * 2);
		const int maximumIterations = static_cast<int>(
			std::min<size_t>(20000, requestedIterations));

		for (; solveIterations < maximumIterations; ++solveIterations) {
			const Eigen::VectorXd gradient =
				system * acceleratedDistances - rhs;
			// Each accepted iterate minimizes the energy while remaining feasible.
			Eigen::VectorXd nextDistances =
				(acceleratedDistances - inverseLipschitz * gradient)
					.cwiseMax(lowerBounds)
					.cwiseMin(upperBounds);

			if (!nextDistances.allFinite()) {
				return depthCells;
			}

			solveError =
				(nextDistances - solvedDistances).norm() /
				std::max(1.0, solvedDistances.norm());
			if (solveError <= 1e-8) {
				solvedDistances = nextDistances;
				++solveIterations;
				break;
			}

			const double nextMomentum =
				0.5 * (1.0 + std::sqrt(1.0 + 4.0 * momentum * momentum));
			acceleratedDistances =
				nextDistances +
				((momentum - 1.0) / nextMomentum) *
					(nextDistances - solvedDistances);
			solvedDistances = nextDistances;
			momentum = nextMomentum;
		}

		// Report the projected-gradient residual of the constrained optimum.
		const Eigen::VectorXd projectedDistances =
			(solvedDistances - inverseLipschitz *
				(system * solvedDistances - rhs))
				.cwiseMax(lowerBounds)
				.cwiseMin(upperBounds);
		solveError =
			(projectedDistances - solvedDistances).norm() /
			std::max(1.0, solvedDistances.norm());
	}

	std::cout << "  biharmonic sparse solve done. Iterations: "
			  << solveIterations
			  << ", error: " << solveError << "\n";
	std::cout.flush();

	// Write the optimized depths back and update their debug classification.
	for (uint i = 0; i < unknownIds.size(); ++i) {
		const uint cellIdx = unknownIds[i];
		const double distance =
			solvedDistances(static_cast<Eigen::Index>(i));

		if (!std::isfinite(distance)) {
			continue;
		}

		CellData& cell = depthCells[cellIdx];
		const bool wasRedHit =
			cell.hasHit && !cell.hasClampedHit;
		const double originalDistance = cell.distance;
		cell.distance = distance;

		// Preserve the legacy collar blend only in the first, unconstrained pass.
		const double originalDistanceWeight =
			originalDistanceWeights[cellIdx];
		if (!useBoxConstraints && originalDistanceWeight > 0.0) {
			cell.distance =
				originalDistanceWeight * originalDistance +
				(1.0 - originalDistanceWeight) * cell.distance;
		}

		cell.hitPoints = {cell.cellCenter + direction * cell.distance};

		// A white cell is cyan when its outside upper bound is active.
		if (useBoxConstraints && !cell.hasHit &&
			std::isfinite(cells[cellIdx].distance)) {
			const double upperBound =
				cells[cellIdx].distance - 0.01 * maxDistance;
			const double activeTolerance = std::max(
				1e-10,
				10.0 * static_cast<double>(eps) *
					std::max(1.0, std::abs(upperBound)));
			cell.isMovedForward =
				cell.distance >= upperBound - activeTolerance;
		}

		if (wasRedHit) {
			cell.isBiharmonicFilledHit = true;
		}
	}

	std::cout << "  biharmonic done\n";
	std::cout.flush();

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
