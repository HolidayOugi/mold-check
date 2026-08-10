#ifndef VCL_TEST_EXTERNAL_888_MOLD_CHECK_BIHARMONIC_H
#define VCL_TEST_EXTERNAL_888_MOLD_CHECK_BIHARMONIC_H

#include "helper.h"
#include "struct.h"

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

struct BiharmonicCellSelection
{
	std::vector<int> variableIds;
	std::vector<unsigned char> fixedIdsMask;
	std::vector<vcl::uint> variableCellIds;
	std::vector<vcl::uint> fixedCellIds;
	std::vector<double> originalDistanceWeights;
};

struct BiharmonicLinearSystem
{
	Eigen::SparseMatrix<double> system;
	Eigen::VectorXd rhs;
};

struct BiharmonicSolveResult
{
	Eigen::VectorXd distances;
	int iterations = 0;
	double error = std::numeric_limits<double>::infinity();
	bool success = false;
};

struct BiharmonicBounds
{
	Eigen::VectorXd lower;
	Eigen::VectorXd upper;
};

static vcl::uint biharmonicNearestHitCellDistance(
	const std::vector<CellData>& depthCells,
	const GridChoice& grid,
	vcl::uint idx,
	vcl::uint radius)
{
	using namespace vcl;

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
			if (!depthCells[neighborIdx].hasHit) {
				continue;
			}

			const uint rowDistance = static_cast<uint>(
				std::abs(r - static_cast<int>(row)));
			const uint colDistance = static_cast<uint>(
				std::abs(c - static_cast<int>(col)));
			nearestDistance = std::min(
				nearestDistance,
				std::max(rowDistance, colDistance));
		}
	}

	return nearestDistance;
}

static bool biharmonicIsAtLeastDistanceFromWhiteBoundary(
	const std::vector<CellData>& depthCells,
	const GridChoice& grid,
	vcl::uint idx,
	vcl::uint minimumDistance)
{
	using namespace vcl;

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
}

static BiharmonicCellSelection biharmonicSelectWhiteFillCells(
	const std::vector<CellData>& cells,
	const std::vector<CellData>& depthCells,
	const GridChoice& grid)
{
	using namespace vcl;

	BiharmonicCellSelection selection;
	selection.variableIds.assign(depthCells.size(), -1);
	selection.fixedIdsMask.assign(depthCells.size(), false);
	selection.variableCellIds.reserve(depthCells.size());

	for (uint idx = 0; idx < depthCells.size(); ++idx) {
		if (!std::isfinite(depthCells[idx].distance)) {
			continue;
		}

		if (cells[idx].hasHit && !cells[idx].hasClampedHit) {
			selection.fixedIdsMask[idx] = true;
			continue;
		}

		if (!cells[idx].hasHit) {
			selection.variableIds[idx] =
				static_cast<int>(selection.variableCellIds.size());
			selection.variableCellIds.push_back(idx);
		}
	}

	for (uint idx = 0; idx < depthCells.size(); ++idx) {
		if (!selection.fixedIdsMask[idx]) {
			continue;
		}

		bool touchesUnknown = false;
		forEachCrossNeighbor(idx, grid, [&](uint neighborIdx) {
			if (selection.variableIds[neighborIdx] >= 0) {
				touchesUnknown = true;
				return false;
			}
			return true;
		});

		if (touchesUnknown) {
			selection.fixedCellIds.push_back(idx);
		}
	}

	return selection;
}

static void biharmonicCollectFixedAnchors(
	const std::vector<CellData>& depthCells,
	const GridChoice& grid,
	vcl::uint collarRadius,
	BiharmonicCellSelection& selection)
{
	using namespace vcl;

	for (uint idx = 0; idx < depthCells.size(); ++idx) {
		if (selection.variableIds[idx] >= 0 ||
			!std::isfinite(depthCells[idx].distance)) {
			continue;
		}

		const uint row = idx / grid.cols;
		const uint col = idx % grid.cols;
		bool touchesUnknown = false;

		if (collarRadius == 0) {
			const auto checkNeighbor = [&](uint neighborIdx) {
				if (selection.variableIds[neighborIdx] >= 0) {
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
					if (selection.variableIds[neighborIdx] >= 0) {
						touchesUnknown = true;
						break;
					}
				}
			}
		}

		if (touchesUnknown) {
			selection.fixedIdsMask[idx] = true;
			selection.fixedCellIds.push_back(idx);
		}
	}
}

static BiharmonicCellSelection biharmonicSelectHitFillCells(
	const std::vector<CellData>& cells,
	std::vector<CellData>& depthCells,
	const GridChoice& grid,
	const vcl::Point3d& direction,
	double eps,
	vcl::uint collarRadius,
	double maxDistance,
	bool useBoxConstraints,
	const std::vector<unsigned char>* fixedWhiteAnchors)
{
	using namespace vcl;

	BiharmonicCellSelection selection;
	selection.variableIds.assign(depthCells.size(), -1);
	selection.fixedIdsMask.assign(depthCells.size(), false);
	selection.originalDistanceWeights.assign(depthCells.size(), 1.0);
	selection.variableCellIds.reserve(depthCells.size());

	const double depthConstraintMargin =
		std::isfinite(maxDistance) ?
			0.003 * maxDistance :
			std::max<double>(1e-6, 100.0 * eps);

	for (uint idx = 0; idx < depthCells.size(); ++idx) {
		if (!std::isfinite(depthCells[idx].distance)) {
			continue;
		}

		const bool isFixedWhiteAnchor =
			useBoxConstraints &&
			fixedWhiteAnchors != nullptr &&
			(*fixedWhiteAnchors)[idx] &&
			!depthCells[idx].hasHit;

		const uint hitDistance =
			depthCells[idx].hasHit ?
				0 :
				biharmonicNearestHitCellDistance(depthCells, grid, idx, collarRadius);

		const bool isWeightedCollar =
			!depthCells[idx].hasHit &&
			collarRadius > 0 &&
			hitDistance <= collarRadius;

		if (depthCells[idx].hasHit) {
			selection.originalDistanceWeights[idx] = 0.0;
		}
		else if (isWeightedCollar) {
			selection.originalDistanceWeights[idx] =
				static_cast<double>(hitDistance) /
				static_cast<double>(collarRadius + 1);
		}

		const bool isFixedOutsideOrange =
			!useBoxConstraints &&
			depthCells[idx].hasHit &&
			!depthCells[idx].hasClampedHit &&
			std::isfinite(cells[idx].distance) &&
			depthCells[idx].distance < cells[idx].distance &&
			biharmonicIsAtLeastDistanceFromWhiteBoundary(
				depthCells,
				grid,
				idx,
				3);

		if (isFixedOutsideOrange) {
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
			isFixedOutsideOrange ||
			depthCells[idx].hasClampedHit ||
			(!depthCells[idx].hasHit && !isWeightedCollar)) {
			continue;
		}

		selection.variableIds[idx] =
			static_cast<int>(selection.variableCellIds.size());
		selection.variableCellIds.push_back(idx);
	}

	biharmonicCollectFixedAnchors(
		depthCells,
		grid,
		collarRadius,
		selection);

	return selection;
}

static BiharmonicLinearSystem biharmonicBuildSystem(
	const std::vector<CellData>& depthCells,
	const GridChoice& grid,
	const std::vector<int>& variableIds,
	const std::vector<unsigned char>& fixedIdsMask,
	const std::vector<vcl::uint>& variableCellIds,
	const std::vector<vcl::uint>& fixedCellIds,
	double eps)
{
	using namespace vcl;

	std::vector<uint> laplacianRowCellIds = variableCellIds;
	laplacianRowCellIds.insert(
		laplacianRowCellIds.end(),
		fixedCellIds.begin(),
		fixedCellIds.end());

	const Eigen::Index rowCount =
		static_cast<Eigen::Index>(laplacianRowCellIds.size());
	const Eigen::Index variableCount =
		static_cast<Eigen::Index>(variableCellIds.size());

	std::vector<Eigen::Triplet<double>> laplacianTriplets;
	laplacianTriplets.reserve(laplacianRowCellIds.size() * 5);
	Eigen::VectorXd fixedRhs = Eigen::VectorXd::Zero(rowCount);

	for (uint rowIdx = 0; rowIdx < laplacianRowCellIds.size(); ++rowIdx) {
		const uint cellIdx = laplacianRowCellIds[rowIdx];
		const uint cellRow = cellIdx / grid.cols;
		const uint cellCol = cellIdx % grid.cols;
		uint usedNeighborCount = 0;

		const auto addNeighbor = [&](uint neighborIdx) {
			if (variableIds[neighborIdx] >= 0) {
				laplacianTriplets.emplace_back(
					static_cast<int>(rowIdx),
					variableIds[neighborIdx],
					1.0);
				++usedNeighborCount;
				return;
			}

			if (fixedIdsMask[neighborIdx]) {
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
			if (variableIds[cellIdx] >= 0) {
				laplacianTriplets.emplace_back(
					static_cast<int>(rowIdx),
					variableIds[cellIdx],
					1.0);
				fixedRhs(static_cast<Eigen::Index>(rowIdx)) =
					depthCells[cellIdx].distance;
			}
			continue;
		}

		if (variableIds[cellIdx] >= 0) {
			laplacianTriplets.emplace_back(
				static_cast<int>(rowIdx),
				variableIds[cellIdx],
				-static_cast<double>(usedNeighborCount));
		}
		else {
			fixedRhs(static_cast<Eigen::Index>(rowIdx)) +=
				static_cast<double>(usedNeighborCount) *
				depthCells[cellIdx].distance;
		}
	}

	Eigen::SparseMatrix<double> laplacian(rowCount, variableCount);
	laplacian.setFromTriplets(
		laplacianTriplets.begin(),
		laplacianTriplets.end());

	BiharmonicLinearSystem linearSystem;
	linearSystem.system = laplacian.transpose() * laplacian;
	linearSystem.rhs = laplacian.transpose() * fixedRhs;

	const double regularization =
		std::max(1e-12, eps * eps);
	for (uint i = 0; i < variableCellIds.size(); ++i) {
		linearSystem.system.coeffRef(
			static_cast<Eigen::Index>(i),
			static_cast<Eigen::Index>(i)) += regularization;
		linearSystem.rhs(static_cast<Eigen::Index>(i)) +=
			regularization * depthCells[variableCellIds[i]].distance;
	}
	linearSystem.system.makeCompressed();

	return linearSystem;
}

static BiharmonicSolveResult biharmonicSolveUnconstrained(
	const BiharmonicLinearSystem& linearSystem,
	size_t variableCount)
{
	BiharmonicSolveResult result;
	result.distances = Eigen::VectorXd(
		static_cast<Eigen::Index>(variableCount));

	Eigen::ConjugateGradient<
		Eigen::SparseMatrix<double>,
		Eigen::Lower | Eigen::Upper> solver;
	solver.setTolerance(1e-8);
	solver.setMaxIterations(
		static_cast<int>(std::max<size_t>(1000, variableCount * 2)));
	solver.compute(linearSystem.system);

	if (solver.info() != Eigen::Success) {
		return result;
	}

	result.distances = solver.solve(linearSystem.rhs);
	result.iterations = solver.iterations();
	result.error = solver.error();
	result.success = solver.info() == Eigen::Success;

	return result;
}

static double biharmonicLipschitzUpperBound(
	const Eigen::SparseMatrix<double>& system)
{
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

	return lipschitzConstant;
}

static BiharmonicSolveResult biharmonicSolveBoxConstrained(
	const BiharmonicLinearSystem& linearSystem,
	const Eigen::VectorXd& initialDistances,
	const BiharmonicBounds& bounds)
{
	BiharmonicSolveResult result;
	result.distances = initialDistances;

	const double lipschitzConstant =
		biharmonicLipschitzUpperBound(linearSystem.system);

	if (!std::isfinite(lipschitzConstant) ||
		lipschitzConstant <= 0.0 ||
		!result.distances.allFinite()) {
		return result;
	}

	Eigen::VectorXd acceleratedDistances = result.distances;
	double momentum = 1.0;
	const double inverseLipschitz = 1.0 / lipschitzConstant;
	const size_t requestedIterations =
		std::max<size_t>(1000, static_cast<size_t>(result.distances.size()) * 2);
	const int maximumIterations = static_cast<int>(
		std::min<size_t>(20000, requestedIterations));

	for (; result.iterations < maximumIterations; ++result.iterations) {
		const Eigen::VectorXd gradient =
			linearSystem.system * acceleratedDistances -
			linearSystem.rhs;
		Eigen::VectorXd nextDistances =
			(acceleratedDistances - inverseLipschitz * gradient)
				.cwiseMax(bounds.lower)
				.cwiseMin(bounds.upper);

		if (!nextDistances.allFinite()) {
			return result;
		}

		result.error =
			(nextDistances - result.distances).norm() /
			std::max(1.0, result.distances.norm());
		if (result.error <= 1e-8) {
			result.distances = nextDistances;
			++result.iterations;
			break;
		}

		const double nextMomentum =
			0.5 * (1.0 + std::sqrt(1.0 + 4.0 * momentum * momentum));
		acceleratedDistances =
			nextDistances +
			((momentum - 1.0) / nextMomentum) *
				(nextDistances - result.distances);
		result.distances = nextDistances;
		momentum = nextMomentum;
	}

	const Eigen::VectorXd projectedDistances =
		(result.distances - inverseLipschitz *
			(linearSystem.system * result.distances - linearSystem.rhs))
			.cwiseMax(bounds.lower)
			.cwiseMin(bounds.upper);
	result.error =
		(projectedDistances - result.distances).norm() /
		std::max(1.0, result.distances.norm());
	result.success = true;

	return result;
}

static BiharmonicBounds biharmonicMakeEmptyBounds(size_t variableCount)
{
	BiharmonicBounds bounds;
	bounds.lower = Eigen::VectorXd::Constant(
		static_cast<Eigen::Index>(variableCount),
		-std::numeric_limits<double>::infinity());
	bounds.upper = Eigen::VectorXd::Constant(
		static_cast<Eigen::Index>(variableCount),
		std::numeric_limits<double>::infinity());

	return bounds;
}

static BiharmonicSolveResult biharmonicSolveWhiteSystem(
	const std::vector<CellData>& cells,
	const std::vector<CellData>& depthCells,
	const std::vector<vcl::uint>& variableCellIds,
	const BiharmonicLinearSystem& linearSystem,
	double maxDistance)
{
	const size_t variableCount = variableCellIds.size();

	if (!std::isfinite(maxDistance)) {
		return biharmonicSolveUnconstrained(linearSystem, variableCount);
	}

	BiharmonicBounds bounds =
		biharmonicMakeEmptyBounds(variableCount);
	Eigen::VectorXd initialDistances(
		static_cast<Eigen::Index>(variableCount));

	for (vcl::uint i = 0; i < variableCellIds.size(); ++i) {
		const vcl::uint cellIdx = variableCellIds[i];
		const double surfaceDistance = cells[cellIdx].distance;
		const double upperBound =
			std::isfinite(surfaceDistance) ?
				surfaceDistance - 0.01 * maxDistance :
				std::numeric_limits<double>::infinity();
		bounds.upper(static_cast<Eigen::Index>(i)) = upperBound;
		initialDistances(static_cast<Eigen::Index>(i)) = std::min(
			depthCells[cellIdx].distance,
			upperBound);
	}

	return biharmonicSolveBoxConstrained(
		linearSystem,
		initialDistances,
		bounds);
}

static BiharmonicBounds biharmonicBuildHitBounds(
	const std::vector<CellData>& cells,
	const std::vector<CellData>& depthCells,
	const std::vector<vcl::uint>& variableCellIds,
	const GridChoice& grid,
	const vcl::Point3d& direction,
	double maxDistance)
{
	BiharmonicBounds bounds =
		biharmonicMakeEmptyBounds(variableCellIds.size());

	const double depthConstraintMargin = 0.003 * maxDistance;

	for (vcl::uint i = 0; i < variableCellIds.size(); ++i) {
		const vcl::uint cellIdx = variableCellIds[i];
		const Eigen::Index variableIdx =
			static_cast<Eigen::Index>(i);
		const CellData& depthCell = depthCells[cellIdx];
		const CellData& surfaceCell = cells[cellIdx];

		if (!depthCell.hasHit) {
			if (std::isfinite(surfaceCell.distance)) {
				bounds.upper(variableIdx) =
					surfaceCell.distance - 0.01 * maxDistance;
			}
		}
		else {
			const bool requiresInteriorDepth =
				biharmonicIsAtLeastDistanceFromWhiteBoundary(
					depthCells,
					grid,
					cellIdx,
					3);
			if (requiresInteriorDepth &&
				std::isfinite(surfaceCell.distance)) {
				bounds.lower(variableIdx) =
					surfaceCell.distance + depthConstraintMargin;
			}

			if (surfaceCell.hitPoints.size() > 1 &&
				std::isfinite(surfaceCell.distance)) {
				const double lastHitDistance =
					(surfaceCell.hitPoints.back() -
					 surfaceCell.cellCenter).dot(direction);
				if (std::isfinite(lastHitDistance) &&
					lastHitDistance > surfaceCell.distance) {
					const double availableMargin = std::min(
						depthConstraintMargin,
						0.5 * (lastHitDistance - surfaceCell.distance));
					bounds.upper(variableIdx) =
						lastHitDistance - availableMargin;
					if (requiresInteriorDepth) {
						bounds.lower(variableIdx) =
							surfaceCell.distance + availableMargin;
					}
				}
			}
		}
	}

	return bounds;
}

static Eigen::VectorXd biharmonicInitialDistances(
	const std::vector<CellData>& depthCells,
	const std::vector<vcl::uint>& variableCellIds,
	const BiharmonicBounds& bounds)
{
	Eigen::VectorXd initialDistances(
		static_cast<Eigen::Index>(variableCellIds.size()));

	for (vcl::uint i = 0; i < variableCellIds.size(); ++i) {
		const Eigen::Index variableIdx =
			static_cast<Eigen::Index>(i);
		double initialDistance = depthCells[variableCellIds[i]].distance;
		initialDistance = std::max(
			initialDistance,
			bounds.lower(variableIdx));
		initialDistance = std::min(
			initialDistance,
			bounds.upper(variableIdx));
		initialDistances(variableIdx) = initialDistance;
	}

	return initialDistances;
}

static void biharmonicApplyWhiteSolution(
	const std::vector<CellData>& cells,
	std::vector<CellData>& depthCells,
	const std::vector<vcl::uint>& variableCellIds,
	const Eigen::VectorXd& solvedDistances,
	const vcl::Point3d& direction,
	double eps,
	double maxDistance)
{
	for (vcl::uint i = 0; i < variableCellIds.size(); ++i) {
		const vcl::uint cellIdx = variableCellIds[i];
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
}

static void biharmonicApplyHitSolution(
	const std::vector<CellData>& cells,
	std::vector<CellData>& depthCells,
	const std::vector<vcl::uint>& variableCellIds,
	const std::vector<double>& originalDistanceWeights,
	const Eigen::VectorXd& solvedDistances,
	const vcl::Point3d& direction,
	double eps,
	double maxDistance,
	bool useBoxConstraints)
{
	for (vcl::uint i = 0; i < variableCellIds.size(); ++i) {
		const vcl::uint cellIdx = variableCellIds[i];
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

		const double originalDistanceWeight =
			originalDistanceWeights[cellIdx];
		if (!useBoxConstraints && originalDistanceWeight > 0.0) {
			cell.distance =
				originalDistanceWeight * originalDistance +
				(1.0 - originalDistanceWeight) * cell.distance;
		}

		cell.hitPoints = {cell.cellCenter + direction * cell.distance};

		if (useBoxConstraints &&
			!cell.hasHit &&
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
}

static std::vector<CellData> biharmonicFillWhiteCellsFromRedCells(
	const std::vector<CellData>& cells,
	std::vector<CellData> depthCells,
	const GridChoice& grid,
	const vcl::Point3d& direction,
	double eps,
	double maxDistance = std::numeric_limits<double>::infinity())
{
	if (depthCells.size() != cells.size() ||
		depthCells.size() != grid.rows * grid.cols) {
		return depthCells;
	}

	BiharmonicCellSelection selection =
		biharmonicSelectWhiteFillCells(cells, depthCells, grid);

	if (selection.variableCellIds.empty() ||
		selection.fixedCellIds.empty()) {
		return depthCells;
	}

	std::cout << "  biharmonic white cells from red cells. Unknown cells: "
			  << selection.variableCellIds.size()
			  << ", fixed anchor cells: " << selection.fixedCellIds.size()
			  << "\n";
	std::cout.flush();

	const BiharmonicLinearSystem linearSystem =
		biharmonicBuildSystem(
			depthCells,
			grid,
			selection.variableIds,
			selection.fixedIdsMask,
			selection.variableCellIds,
			selection.fixedCellIds,
			eps);

	std::cout << "  biharmonic sparse solve start\n";
	std::cout.flush();

	const BiharmonicSolveResult solveResult =
		biharmonicSolveWhiteSystem(
			cells,
			depthCells,
			selection.variableCellIds,
			linearSystem,
			maxDistance);

	if (!solveResult.success) {
		return depthCells;
	}

	std::cout << "  biharmonic sparse solve done. Iterations: "
			  << solveResult.iterations
			  << ", error: " << solveResult.error << "\n";
	std::cout.flush();

	biharmonicApplyWhiteSolution(
		cells,
		depthCells,
		selection.variableCellIds,
		solveResult.distances,
		direction,
		eps,
		maxDistance);

	std::cout << "  biharmonic done\n";
	std::cout.flush();

	return depthCells;
}

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
	if (depthCells.size() != cells.size() ||
		depthCells.size() != grid.rows * grid.cols) {
		return depthCells;
	}

	const bool useBoxConstraints =
		fixedWhiteAnchors != nullptr &&
		fixedWhiteAnchors->size() == depthCells.size() &&
		std::isfinite(maxDistance);

	BiharmonicCellSelection selection =
		biharmonicSelectHitFillCells(
			cells,
			depthCells,
			grid,
			direction,
			eps,
			collarRadius,
			maxDistance,
			useBoxConstraints,
			fixedWhiteAnchors);

	if (selection.variableCellIds.empty() ||
		selection.fixedCellIds.empty()) {
		return depthCells;
	}

	std::cout << "  biharmonic unknown cells: "
			  << selection.variableCellIds.size()
			  << ", fixed anchor cells: " << selection.fixedCellIds.size()
			  << "\n";
	std::cout.flush();

	const BiharmonicLinearSystem linearSystem =
		biharmonicBuildSystem(
			depthCells,
			grid,
			selection.variableIds,
			selection.fixedIdsMask,
			selection.variableCellIds,
			selection.fixedCellIds,
			eps);

	std::cout << "  biharmonic sparse solve start\n";
	std::cout.flush();

	BiharmonicSolveResult solveResult;
	if (!useBoxConstraints) {
		solveResult = biharmonicSolveUnconstrained(
			linearSystem,
			selection.variableCellIds.size());
	}
	else {
		const BiharmonicBounds bounds =
			biharmonicBuildHitBounds(
				cells,
				depthCells,
				selection.variableCellIds,
				grid,
				direction,
				maxDistance);
		const Eigen::VectorXd initialDistances =
			biharmonicInitialDistances(
				depthCells,
				selection.variableCellIds,
				bounds);
		solveResult = biharmonicSolveBoxConstrained(
			linearSystem,
			initialDistances,
			bounds);
	}

	if (!solveResult.success) {
		return depthCells;
	}

	std::cout << "  biharmonic sparse solve done. Iterations: "
			  << solveResult.iterations
			  << ", error: " << solveResult.error << "\n";
	std::cout.flush();

	biharmonicApplyHitSolution(
		cells,
		depthCells,
		selection.variableCellIds,
		selection.originalDistanceWeights,
		solveResult.distances,
		direction,
		eps,
		maxDistance,
		useBoxConstraints);

	std::cout << "  biharmonic done\n";
	std::cout.flush();

	return depthCells;
}
#endif
