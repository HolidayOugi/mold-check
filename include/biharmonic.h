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

static void biharmonicSetCellDistance(
	CellData& cell,
	double distance,
	const vcl::Point3d& direction)
{
	cell.distance = distance;
	const vcl::Point3d hitPoint =
		cell.cellCenter + direction * cell.distance;
	if (cell.hitPoints.empty()) {
		cell.hitPoints.push_back(hitPoint);
	}
	else {
		cell.hitPoints[0] = hitPoint;
	}
}

// Recompute the current inside/outside state from the first mesh hit only.
static void updateDepthCellInsideFlags(
	const std::vector<CellData>& cells,
	std::vector<CellData>& depthCells,
	double eps)
{
	if (cells.size() != depthCells.size()) {
		return;
	}

	const double tolerance = std::max(1e-10, 10.0 * std::abs(eps));
	for (vcl::uint idx = 0; idx < depthCells.size(); ++idx) {
		const CellData& surfaceCell = cells[idx];
		CellData& depthCell = depthCells[idx];

		depthCell.isInside = false;
		if (!surfaceCell.hasHit ||
			!std::isfinite(surfaceCell.distance) ||
			!std::isfinite(depthCell.distance)) {
			continue;
		}

		depthCell.isInside =
			depthCell.distance >= surfaceCell.distance - tolerance;
	}
}

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

// Select white cells as unknowns and original hits as fixed anchors.
static BiharmonicCellSelection biharmonicSelectWhiteFillCells(
	const std::vector<CellData>& cells,
	const std::vector<CellData>& depthCells,
	const GridChoice& grid)
{
	using namespace vcl;

	// Start with no variables or anchors selected.
	BiharmonicCellSelection selection;
	selection.variableIds.assign(depthCells.size(), -1);
	selection.fixedIdsMask.assign(depthCells.size(), false);
	selection.variableCellIds.reserve(depthCells.size());

	for (uint idx = 0; idx < depthCells.size(); ++idx) {
		if (depthCells[idx].isDiscarded ||
			!std::isfinite(depthCells[idx].distance)) {
			continue;
		}

		// Use non-clamped original hits as fixed cells.
		if (cells[idx].hasHit && !cells[idx].hasClampedHit) {
			selection.fixedIdsMask[idx] = true;
			continue;
		}

		// Map each white cell to one solver variable.
		if (!cells[idx].hasHit) {
			selection.variableIds[idx] =
				static_cast<int>(selection.variableCellIds.size());
			selection.variableCellIds.push_back(idx);
		}
	}

	// Keep only fixed cells touching an unknown as explicit anchors.
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
		if (depthCells[idx].isDiscarded ||
			selection.variableIds[idx] >= 0 ||
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

// Select orange/hit cells and the movable white/cyan collar.
static BiharmonicCellSelection biharmonicSelectHitFillCells(
	const std::vector<CellData>& cells,
	std::vector<CellData>& depthCells,
	const GridChoice& grid,
	const vcl::Point3d& direction,
	double eps,
	double maxDistance,
	bool useBoxConstraints,
	vcl::uint collarRadius)
{
	using namespace vcl;

	// Start with no variables or anchors selected.
	BiharmonicCellSelection selection;
	selection.variableIds.assign(depthCells.size(), -1);
	selection.fixedIdsMask.assign(depthCells.size(), false);
	selection.originalDistanceWeights.assign(depthCells.size(), 0.0);
	selection.variableCellIds.reserve(depthCells.size());

	// Small safety margin used when an orange point must stay inside the mesh.
	const double depthConstraintMargin =
		std::isfinite(maxDistance) ?
			0.003 * maxDistance :
			std::max<double>(1e-6, 100.0 * eps);

	// Refresh inside/outside state before deciding which hit cells can move.
	updateDepthCellInsideFlags(cells, depthCells, eps);

	for (uint idx = 0; idx < depthCells.size(); ++idx) {
		if (depthCells[idx].isDiscarded ||
			!std::isfinite(depthCells[idx].distance)) {
			continue;
		}

		const uint hitDistance =
			depthCells[idx].hasHit ?
				0 :
				biharmonicNearestHitCellDistance(
					depthCells,
					grid,
					idx,
					collarRadius);

		// Select non-hit cells inside the collar when box bounds are active.
		const bool isCollarCell =
			useBoxConstraints &&
			!depthCells[idx].hasHit &&
			collarRadius > 0 &&
			hitDistance <= collarRadius;

		// Build a smooth blend weight: near-hit cells follow the solution,
		// while cells near the outer collar keep more of their original depth.
		if (isCollarCell) {
			const double normalizedDistance =
				static_cast<double>(hitDistance) /
				static_cast<double>(collarRadius + 1);
			selection.originalDistanceWeights[idx] =
				normalizedDistance * normalizedDistance *
				(3.0 - 2.0 * normalizedDistance);
		}

		// Deep orange points that accidentally moved outside are fixed back inside.
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
			// Move it just after the first hit, but avoid crossing the last hit if present.
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

			biharmonicSetCellDistance(
				constrainedCell,
				constrainedDistance,
				direction);
			constrainedCell.isBiharmonicFilledHit = true;
		}

		// In the second pass, both white and cyan collar cells join the solve.
		if (isFixedOutsideOrange ||
			depthCells[idx].hasClampedHit ||
			(!depthCells[idx].hasHit && !isCollarCell)) {
			continue;
		}

		// Everything left is a variable solved by the hit biharmonic pass.
		selection.variableIds[idx] =
			static_cast<int>(selection.variableCellIds.size());
		selection.variableCellIds.push_back(idx);
	}

	// Some points may have been corrected, so refresh before collecting anchors.
	updateDepthCellInsideFlags(cells, depthCells, eps);

	// Add the fixed band outside the movable collar.
	biharmonicCollectFixedAnchors(
		depthCells,
		grid,
		collarRadius,
		selection);

	return selection;
}

// Build the shared biharmonic least-squares system from variables and anchors.
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

	// Build Laplacian rows for unknowns and boundary anchors.
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

		// Add unknown neighbors to the matrix and fixed ones to the RHS.
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
			// Keep isolated variables at their current distance.
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

	// Form the normal equations for the Laplacian least-squares problem.
	BiharmonicLinearSystem linearSystem;
	linearSystem.system = laplacian.transpose() * laplacian;
	linearSystem.rhs = laplacian.transpose() * fixedRhs;

	// Regularize the system around the current distances.
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

// Fast path: solve the unconstrained biharmonic system with Eigen CG.
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

// Estimate a safe step size for the projected-gradient solver.
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

// Bound path: project each step into the hard distance bounds.
static BiharmonicSolveResult biharmonicSolveBoxConstrained(
	const BiharmonicLinearSystem& linearSystem,
	const Eigen::VectorXd& initialDistances,
	const BiharmonicBounds& bounds)
{
	BiharmonicSolveResult result;
	result.distances = initialDistances;

	// The Lipschitz value keeps each gradient step stable.
	const double lipschitzConstant =
		biharmonicLipschitzUpperBound(linearSystem.system);

	if (!std::isfinite(lipschitzConstant) ||
		lipschitzConstant <= 0.0 ||
		!result.distances.allFinite()) {
		return result;
	}

	// FISTA-style acceleration starts from the supplied initial guess.
	Eigen::VectorXd acceleratedDistances = result.distances;
	double momentum = 1.0;
	const double inverseLipschitz = 1.0 / lipschitzConstant;
	// Match the normal CG budget scale, then cap it to avoid runaway solves.
	const size_t requestedIterations =
		std::max<size_t>(1000, static_cast<size_t>(result.distances.size()) * 2);
	const int maximumIterations = static_cast<int>(
		std::min<size_t>(30000, requestedIterations));

	for (; result.iterations < maximumIterations; ++result.iterations) {
		// Take one descent step and immediately project it into the hard bounds.
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

		// Stop when the projected solution stops changing enough.
		result.error =
			(nextDistances - result.distances).norm() /
			std::max(1.0, result.distances.norm());
		if (result.error <= 1e-8) {
			result.distances = nextDistances;
			++result.iterations;
			break;
		}

		// Update the acceleration term for the next projected step.
		const double nextMomentum =
			0.5 * (1.0 + std::sqrt(1.0 + 4.0 * momentum * momentum));
		acceleratedDistances =
			nextDistances +
			((momentum - 1.0) / nextMomentum) *
				(nextDistances - result.distances);
		result.distances = nextDistances;
		momentum = nextMomentum;
	}

	// Report the final projected-step error, not only the last iteration change.
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

// Tunable height offset parameters, expressed as fractions of maxDistance.
static constexpr double BiharmonicWhiteHeightMaxFraction = 0.15;
static constexpr double BiharmonicWhiteHeightGrowthPerCellFraction = 0.005;

static vcl::uint biharmonicInvalidBoundaryDistance()
{
	return std::numeric_limits<vcl::uint>::max();
}

static bool biharmonicIsWhiteForwardCapCandidate(
	const std::vector<CellData>& cells,
	const std::vector<CellData>& depthCells,
	vcl::uint idx)
{
	return idx < cells.size() &&
		idx < depthCells.size() &&
		!cells[idx].hasHit &&
		cells[idx].isReduced &&
		!depthCells[idx].hasHit &&
		std::isfinite(cells[idx].distance) &&
		std::isfinite(depthCells[idx].distance);
}

static std::vector<vcl::uint> biharmonicWhiteBoundaryDistances(
	const std::vector<CellData>& cells,
	const std::vector<CellData>& depthCells,
	const GridChoice& grid)
{
	using namespace vcl;

	const uint invalidDistance = biharmonicInvalidBoundaryDistance();
	std::vector<uint> distances(depthCells.size(), invalidDistance);
	std::vector<uint> queue;
	queue.reserve(depthCells.size());

	const auto seedBoundaryWhites = [&](bool orangeOnly) {
		for (uint idx = 0; idx < depthCells.size(); ++idx) {
			if (distances[idx] != invalidDistance ||
				!biharmonicIsWhiteForwardCapCandidate(
					cells,
					depthCells,
					idx)) {
				continue;
			}

			bool touchesBoundaryHit = false;
			forEachCrossNeighbor(idx, grid, [&](uint neighborIdx) {
				const bool isBoundaryHit =
					neighborIdx < depthCells.size() &&
					depthCells[neighborIdx].hasHit &&
					(!orangeOnly || depthCells[neighborIdx].isBiharmonicFilledHit);
				if (isBoundaryHit) {
					touchesBoundaryHit = true;
					return false;
				}
				return true;
			});

			if (touchesBoundaryHit) {
				// The first white-cell ring next to the orange boundary has d = 1.
				distances[idx] = 1;
				queue.push_back(idx);
			}
		}
	};

	seedBoundaryWhites(true);
	if (queue.empty()) {
		seedBoundaryWhites(false);
	}

	for (size_t readIdx = 0; readIdx < queue.size(); ++readIdx) {
		const uint idx = queue[readIdx];
		forEachCrossNeighbor(idx, grid, [&](uint neighborIdx) {
			if (distances[neighborIdx] == invalidDistance &&
				biharmonicIsWhiteForwardCapCandidate(
					cells,
					depthCells,
					neighborIdx)) {
				distances[neighborIdx] = distances[idx] + 1;
				queue.push_back(neighborIdx);
			}
			return true;
		});
	}

	return distances;
}


// Compute upperBound = originalDistance - min(hMax, d * heightGrowthPerCell).
static double biharmonicWhiteForwardCapDistance(
	const std::vector<CellData>& cells,
	vcl::uint cellIdx,
	double maxDistance,
	const std::vector<vcl::uint>& whiteBoundaryDistances)
{
	if (cellIdx >= cells.size() ||
		cellIdx >= whiteBoundaryDistances.size() ||
		!std::isfinite(cells[cellIdx].distance) ||
		!std::isfinite(maxDistance)) {
		return std::numeric_limits<double>::infinity();
	}

	const vcl::uint invalidDistance = biharmonicInvalidBoundaryDistance();
	const vcl::uint boundaryDistance = whiteBoundaryDistances[cellIdx];
	if (boundaryDistance == invalidDistance) {
		return std::numeric_limits<double>::infinity();
	}

	const double maximumHeight =
		BiharmonicWhiteHeightMaxFraction * maxDistance;
	const double heightGrowthPerCell =
		BiharmonicWhiteHeightGrowthPerCellFraction * maxDistance;
	const double height = std::min(
		maximumHeight,
		static_cast<double>(boundaryDistance) * heightGrowthPerCell);
	return cells[cellIdx].distance - height;
}
// Constrain white cells using the distance-dependent height offset.
static BiharmonicSolveResult biharmonicSolveWhiteSystem(
	const std::vector<CellData>& cells,
	const std::vector<CellData>& depthCells,
	const std::vector<vcl::uint>& variableCellIds,
	const GridChoice& grid,
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
	const std::vector<vcl::uint> whiteBoundaryDistances =
		biharmonicWhiteBoundaryDistances(cells, depthCells, grid);

	for (vcl::uint i = 0; i < variableCellIds.size(); ++i) {
		const vcl::uint cellIdx = variableCellIds[i];
		const Eigen::Index variableIdx =
			static_cast<Eigen::Index>(i);
		initialDistances(variableIdx) = depthCells[cellIdx].distance;

		if (biharmonicIsWhiteForwardCapCandidate(
				cells,
				depthCells,
				cellIdx)) {
			bounds.upper(variableIdx) =
				biharmonicWhiteForwardCapDistance(
					cells,
					cellIdx,
					maxDistance,
					whiteBoundaryDistances);
		}
	}

	return biharmonicSolveBoxConstrained(
		linearSystem,
		initialDistances,
		bounds);
}
// Build hard bounds that keep inside orange points safely inside the mesh.
static BiharmonicBounds biharmonicBuildHitBounds(
	const std::vector<CellData>& cells,
	const std::vector<CellData>& depthCells,
	const std::vector<vcl::uint>& variableCellIds,
	const GridChoice& grid,
	const vcl::Point3d& direction,
	double maxDistance,
	const std::vector<unsigned char>* cyanCells)
{
	BiharmonicBounds bounds =
		biharmonicMakeEmptyBounds(variableCellIds.size());
	const std::vector<vcl::uint> whiteBoundaryDistances =
		biharmonicWhiteBoundaryDistances(cells, depthCells, grid);

	// Orange points that are clearly inside get this safety margin after the first hit.
	const double depthConstraintMargin = 0.003 * maxDistance;

	for (vcl::uint i = 0; i < variableCellIds.size(); ++i) {
		const vcl::uint cellIdx = variableCellIds[i];
		const Eigen::Index variableIdx =
			static_cast<Eigen::Index>(i);
		const CellData& depthCell = depthCells[cellIdx];
		const CellData& surfaceCell = cells[cellIdx];

		if (!depthCell.hasHit) {
			// Reuse the same distance-minus-height upper bound in the hit pass.
			if (biharmonicIsWhiteForwardCapCandidate(
					cells,
					depthCells,
					cellIdx)) {
				bounds.upper(variableIdx) =
					biharmonicWhiteForwardCapDistance(
						cells,
						cellIdx,
						maxDistance,
						whiteBoundaryDistances);
			}

			// Cyan cells may participate in the collar, but never increase
			// their distance relative to the value entering this hit pass.
			if (cyanCells != nullptr &&
				cyanCells->size() == depthCells.size() &&
				(*cyanCells)[cellIdx]) {
				bounds.upper(variableIdx) = std::min(
					bounds.upper(variableIdx),
					depthCell.distance);
			}
		}
		else {
			// Orange cells at least three cells from white boundary must stay inside.
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

			// If the ray exits the mesh later, keep the value before that last hit.
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

// Start the constrained solve from current values already projected into bounds.
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
		// Clamp the starting value so the solver begins from a valid state.
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

// Write solved distances back and mark whites that reached their computed
// originalDistance-minus-height upper bound as cyan.
static void biharmonicApplyWhiteSolution(
	const std::vector<CellData>& cells,
	std::vector<CellData>& depthCells,
	const std::vector<vcl::uint>& variableCellIds,
	const GridChoice& grid,
	const Eigen::VectorXd& solvedDistances,
	const vcl::Point3d& direction,
	double eps,
	double maxDistance)
{
	const std::vector<vcl::uint> whiteBoundaryDistances =
		std::isfinite(maxDistance) ?
			biharmonicWhiteBoundaryDistances(cells, depthCells, grid) :
			std::vector<vcl::uint>();

	//set distance to each cell
	for (vcl::uint i = 0; i < variableCellIds.size(); ++i) {
		const vcl::uint cellIdx = variableCellIds[i];
		const double distance =
			solvedDistances(static_cast<Eigen::Index>(i));

		if (!std::isfinite(distance)) {
			continue;
		}

		CellData& cell = depthCells[cellIdx];
		biharmonicSetCellDistance(cell, distance, direction);
		
		// Mark cells touching their distance-minus-height cap as cyan.
		if (std::isfinite(maxDistance) &&
			std::isfinite(cells[cellIdx].distance) &&
			biharmonicIsWhiteForwardCapCandidate(
				cells,
				depthCells,
				cellIdx)) {
			const double upperBound =
				biharmonicWhiteForwardCapDistance(
					cells,
					cellIdx,
					maxDistance,
					whiteBoundaryDistances);
			const double activeTolerance = std::max(
				1e-10,
				10.0 * static_cast<double>(eps) *
					std::max(1.0, std::abs(upperBound)));
			cell.isMovedForward =
				distance >= upperBound - activeTolerance;
		}
	}

	updateDepthCellInsideFlags(cells, depthCells, eps);
}

// Write solved hit distances back and refresh color/state flags.
static void biharmonicApplyHitSolution(
	const std::vector<CellData>& cells,
	std::vector<CellData>& depthCells,
	const std::vector<vcl::uint>& variableCellIds,
	const std::vector<double>& originalDistanceWeights,
	const std::vector<unsigned char>* cyanCells,
	const GridChoice& grid,
	const Eigen::VectorXd& solvedDistances,
	const vcl::Point3d& direction,
	double eps,
	double maxDistance,
	bool useBoxConstraints)
{
	// Precompute how far each cell lies from the white boundary.
	const std::vector<vcl::uint> whiteBoundaryDistances =
		useBoxConstraints ?
			biharmonicWhiteBoundaryDistances(cells, depthCells, grid) :
			std::vector<vcl::uint>();

	// Write each valid solver result back to its grid cell.
	for (vcl::uint i = 0; i < variableCellIds.size(); ++i) {
		const vcl::uint cellIdx = variableCellIds[i];
		const double distance =
			solvedDistances(static_cast<Eigen::Index>(i));

		if (!std::isfinite(distance)) {
			continue;
		}

		// Remember the incoming state before replacing the cell depth.
		CellData& cell = depthCells[cellIdx];
		const bool wasCyan =
			cyanCells != nullptr &&
			cyanCells->size() == depthCells.size() &&
			(*cyanCells)[cellIdx];
		const bool wasRedHit =
			cell.hasHit && !cell.hasClampedHit;
		const double originalDistance = cell.distance;
		biharmonicSetCellDistance(cell, distance, direction);

		// Reduce motion toward the outer edge of the movable white collar.
		if (useBoxConstraints &&
			!cell.hasHit &&
			cellIdx < originalDistanceWeights.size()) {
			const double originalWeight = std::clamp(
				originalDistanceWeights[cellIdx],
				0.0,
				1.0);
			biharmonicSetCellDistance(
				cell,
				originalWeight * originalDistance +
					(1.0 - originalWeight) * cell.distance,
				direction);
		}

		// Preserve existing cyan cells and mark cells still touching their
		// distance-minus-height cap.
		if (useBoxConstraints &&
			!cell.hasHit &&
			std::isfinite(cells[cellIdx].distance) &&
			biharmonicIsWhiteForwardCapCandidate(
				cells,
				depthCells,
				cellIdx)) {
			const double upperBound =
				biharmonicWhiteForwardCapDistance(
					cells,
					cellIdx,
					maxDistance,
					whiteBoundaryDistances);
			const double activeTolerance = std::max(
				1e-10,
				10.0 * static_cast<double>(eps) *
					std::max(1.0, std::abs(upperBound)));
			cell.isMovedForward =
				wasCyan ||
				cell.distance >= upperBound - activeTolerance;
		}
		else if (wasCyan) {
			cell.isMovedForward = true;
		}

		// Mark original red hits updated by this pass.
		if (wasRedHit) {
			cell.isBiharmonicFilledHit = true;
		}
	}

	// Refresh inside/outside flags after all depths have changed.
	updateDepthCellInsideFlags(cells, depthCells, eps);
}

static std::vector<CellData> biharmonicFillWhiteCells(
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

	// Keep inside flags coherent before selecting unknowns and anchors.
	updateDepthCellInsideFlags(cells, depthCells, eps);

	// Pick white variables and fixed hit anchors for this pass.
	BiharmonicCellSelection selection =
		biharmonicSelectWhiteFillCells(cells, depthCells, grid);

	if (selection.variableCellIds.empty() ||
		selection.fixedCellIds.empty()) {
		return depthCells;
	}

	std::cout << "  biharmonic white cells. Unknown cells: "
			  << selection.variableCellIds.size()
			  << ", fixed anchor cells: " << selection.fixedCellIds.size()
			  << "\n";
	std::cout.flush();

	// Build the same biharmonic system used by both passes.
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
			grid,
			linearSystem,
			maxDistance);

	if (!solveResult.success) {
		return depthCells;
	}

	std::cout << "  biharmonic sparse solve done. Iterations: "
			  << solveResult.iterations
			  << ", error: " << solveResult.error << "\n";
	std::cout.flush();

	// Store the solution and update debug/state flags.
	biharmonicApplyWhiteSolution(
		cells,
		depthCells,
		selection.variableCellIds,
		grid,
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
	vcl::uint collarRadius,
	double maxDistance = std::numeric_limits<double>::infinity(),
	const std::vector<unsigned char>* cyanCells = nullptr)
{
	if (depthCells.size() != cells.size() ||
		depthCells.size() != grid.rows * grid.cols) {
		return depthCells;
	}

	// Recompute inside state because orange cells may now be outside or inside.
	updateDepthCellInsideFlags(cells, depthCells, eps);

	// The second hit pass enables the collar and cyan one-sided bounds.
	const bool useBoxConstraints =
		cyanCells != nullptr &&
		cyanCells->size() == depthCells.size() &&
		std::isfinite(maxDistance);
	const vcl::uint activeCollarRadius =
		useBoxConstraints ? collarRadius : 0;

	// Pick movable hit cells and fixed anchors.
	BiharmonicCellSelection selection =
		biharmonicSelectHitFillCells(
			cells,
			depthCells,
			grid,
			direction,
			eps,
			maxDistance,
			useBoxConstraints,
			activeCollarRadius);

	if (selection.variableCellIds.empty() ||
		selection.fixedCellIds.empty()) {
		return depthCells;
	}

	std::cout << "  biharmonic unknown cells: "
			  << selection.variableCellIds.size()
			  << ", fixed anchor cells: " << selection.fixedCellIds.size()
			  << "\n";
	std::cout.flush();

	// Build the same biharmonic system used by both passes.
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

	// Use plain CG unless hard bounds are required by the second pass.
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
				maxDistance,
				cyanCells);
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

	// Store hit/orange results and refresh color/state flags.
	biharmonicApplyHitSolution(
		cells,
		depthCells,
		selection.variableCellIds,
		selection.originalDistanceWeights,
		cyanCells,
		grid,
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
