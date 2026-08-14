#ifndef VCL_TEST_EXTERNAL_888_MOLD_CHECK_REDUCE_H
#define VCL_TEST_EXTERNAL_888_MOLD_CHECK_REDUCE_H

#include "debug_output.h"
#include "helper.h"
#include "struct.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include <vclib/base/parallel.h>
#include <vclib/meshes.h>

static void erodeHitMaskOnce(
	std::vector<CellData>& cells,
	const GridChoice& grid)
{
	using namespace vcl;

	std::vector<uint> allCells(cells.size());
	std::iota(allCells.begin(), allCells.end(), 0);
	std::vector<char> nextHasHit(cells.size(), false);

	parallelFor(allCells, [&](uint idx) {
		if (!cells[idx].hasHit) {
			return;
		}

		bool keep = true;
		forEachCrossNeighbor(idx, grid, [&](uint neighborIdx) {
			if (!cells[neighborIdx].hasHit) {
				keep = false;
				return false;
			}
			return true;
		});

		nextHasHit[idx] = keep;
	});

	parallelFor(allCells, [&](uint idx) {
		if (cells[idx].hasHit && !nextHasHit[idx]) {
			cells[idx].isReduced = true;
		}
		cells[idx].hasHit = nextHasHit[idx];
	});
}

static void dilateHitMaskOnce(
	std::vector<CellData>& cells,
	const GridChoice& grid)
{
	using namespace vcl;

	std::vector<uint> allCells(cells.size());
	std::iota(allCells.begin(), allCells.end(), 0);
	std::vector<char> nextHasHit(cells.size(), false);

	parallelFor(allCells, [&](uint idx) {
		bool add = cells[idx].hasHit;
		forEachCrossNeighbor(idx, grid, [&](uint neighborIdx) {
			if (cells[neighborIdx].hasHit) {
				add = true;
				return false;
			}
			return true;
		});

		nextHasHit[idx] = add;
	});

	parallelFor(allCells, [&](uint idx) {
		if (cells[idx].hasHit && !nextHasHit[idx]) {
			cells[idx].isReduced = true;
		}
		cells[idx].hasHit = nextHasHit[idx];
	});
}

static std::vector<std::vector<vcl::uint>> removeAngleJumpPoints(
	const std::vector<CellData>& cells,
	const GridChoice& grid,
	const vcl::Point3d& direction,
	double angleThresholdDegrees,
	float eps)
{
	using namespace vcl;

	Point3d normalizedDirection = direction;
	double coneCosThreshold = -std::numeric_limits<double>::infinity();
	bool useAngleThreshold = false;

	if (angleThresholdDegrees > 0.0 && normalizedDirection.norm() > eps) {
		normalizedDirection.normalize();
		coneCosThreshold =
			std::cos(angleThresholdDegrees * M_PI / 180.0);
		useAngleThreshold = true;
	}

	std::vector<std::vector<uint>> connectedNeighbors(cells.size());
	for (uint idx = 0; idx < cells.size(); ++idx) {
		if (!cells[idx].hasHit || cells[idx].hitPoints.empty()) {
			continue;
		}

		forEachCrossNeighbor(idx, grid, [&](uint neighborIdx) {
			if (neighborIdx <= idx ||
				!cells[neighborIdx].hasHit ||
				cells[neighborIdx].hitPoints.empty()) {
				return true;
			}

			const Point3d& point = cells[idx].hitPoints[0];
			const Point3d& neighborPoint =
				cells[neighborIdx].hitPoints[0];

			const bool isAngleJump =
				useAngleThreshold &&
				(isWithinPlaneAngle(
					 point,
					 neighborPoint,
					 normalizedDirection,
					 coneCosThreshold,
					 eps));

			if (!isAngleJump) {
				connectedNeighbors[idx].push_back(neighborIdx);
				connectedNeighbors[neighborIdx].push_back(idx);
			}

			return true;
		});
	}

	return connectedNeighbors;
}

static std::vector<CellData> keepLargestHitComponent(
	const std::vector<CellData>& cells,
	const std::vector<std::vector<vcl::uint>>& connectedNeighbors)
{
	using namespace vcl;

	if (cells.size() != connectedNeighbors.size()) {
		return cells;
	}

	std::vector<bool> visited(cells.size(), false);
	std::vector<uint> largestComponent;

	for (uint start = 0; start < cells.size(); ++start) {
		if (visited[start] || !cells[start].hasHit) {
			continue;
		}

		std::vector<uint> component;
		std::vector<uint> stack;

		visited[start] = true;
		stack.push_back(start);

		while (!stack.empty()) {
			const uint idx = stack.back();
			stack.pop_back();
			component.push_back(idx);

			for (uint neighborIdx : connectedNeighbors[idx]) {
				if (!visited[neighborIdx] && cells[neighborIdx].hasHit) {
					visited[neighborIdx] = true;
					stack.push_back(neighborIdx);
				}
			}
		}

		if (component.size() > largestComponent.size()) {
			largestComponent = std::move(component);
		}
	}

	if (largestComponent.empty()) {
		return cells;
	}

	std::vector<char> keep(cells.size(), false);
	for (uint idx : largestComponent) {
		keep[idx] = true;
	}

	std::vector<CellData> result = cells;

	for (uint idx = 0; idx < result.size(); ++idx) {
		if (result[idx].hasHit && !keep[idx]) {
			result[idx].hasHit = false;
			result[idx].isReduced = true;
		}
	}

	return result;
}

static bool isBoundaryHitCell(
	const std::vector<CellData>& cells,
	const GridChoice& grid,
	vcl::uint idx)
{
	using namespace vcl;

	if (!cells[idx].hasHit) {
		return false;
	}

	const uint row = idx / grid.cols;
	const uint col = idx % grid.cols;

	if (row == 0 || col == 0 ||
		row + 1 == grid.rows || col + 1 == grid.cols) {
		return true;
	}

	bool hasEmptyNeighbor = false;
	forEachCrossNeighbor(idx, grid, [&](uint neighborIdx) {
		if (!cells[neighborIdx].hasHit) {
			hasEmptyNeighbor = true;
			return false;
		}
		return true;
	});

	return hasEmptyNeighbor;
}

static std::vector<vcl::uint> boundaryHitCells(
	const std::vector<CellData>& cells,
	const GridChoice& grid)
{
	using namespace vcl;

	std::vector<uint> boundary;
	for (uint idx = 0; idx < cells.size(); ++idx) {
		if (isBoundaryHitCell(cells, grid, idx)) {
			boundary.push_back(idx);
		}
	}

	return boundary;
}

static std::vector<vcl::uint> rasterizedGridSegment(
	vcl::uint a,
	vcl::uint b,
	const GridChoice& grid)
{
	using namespace vcl;

	int row0 = static_cast<int>(a / grid.cols);
	int col0 = static_cast<int>(a % grid.cols);
	const int row1 = static_cast<int>(b / grid.cols);
	const int col1 = static_cast<int>(b % grid.cols);

	const int dCol = std::abs(col1 - col0);
	const int dRow = -std::abs(row1 - row0);
	const int stepCol = (col0 < col1) ? 1 : -1;
	const int stepRow = (row0 < row1) ? 1 : -1;
	int err = dCol + dRow;

	std::vector<uint> line;
	line.reserve(
		static_cast<size_t>(
			std::max(std::abs(col1 - col0), std::abs(row1 - row0)) + 1));
	while (true) {
		if (row0 >= 0 && col0 >= 0 &&
			row0 < static_cast<int>(grid.rows) &&
			col0 < static_cast<int>(grid.cols)) {
			line.push_back(
				static_cast<uint>(row0) * grid.cols +
				static_cast<uint>(col0));
		}

		if (row0 == row1 && col0 == col1) {
			break;
		}

		const int err2 = 2 * err;
		if (err2 >= dRow) {
			err += dRow;
			col0 += stepCol;
		}
		if (err2 <= dCol) {
			err += dCol;
			row0 += stepRow;
		}
	}

	return line;
}

static bool lineInsideHitMask(
	const std::vector<CellData>& cells,
	const std::vector<vcl::uint>& line)
{
	if (line.size() < 2) {
		return false;
	}

	for (vcl::uint idx : line) {
		if (!cells[idx].hasHit) {
			return false;
		}
	}

	return true;
}

static double gridCellCenterDistance(
	vcl::uint a,
	vcl::uint b,
	const GridChoice& grid)
{
	const double rowA = static_cast<double>(a / grid.cols);
	const double colA = static_cast<double>(a % grid.cols);
	const double rowB = static_cast<double>(b / grid.cols);
	const double colB = static_cast<double>(b % grid.cols);

	const double du = (colA - colB) * grid.sideU;
	const double dv = (rowA - rowB) * grid.sideV;
	return std::sqrt(du * du + dv * dv);
}

static HitCellShapeData largestUnblockedHitComponentShape(
	const std::vector<CellData>& cells,
	const GridChoice& grid,
	const std::vector<char>& blocked,
	std::vector<char>* keepMask = nullptr)
{
	using namespace vcl;

	HitCellShapeData result;

	if (blocked.size() != cells.size()) {
		return result;
	}

	std::vector<char> visited(cells.size(), false);
	std::vector<uint> stack;
	std::vector<uint> component;
	std::vector<uint> largestComponent;

	for (uint start = 0; start < cells.size(); ++start) {
		if (visited[start] || blocked[start] || !cells[start].hasHit) {
			continue;
		}

		component.clear();
		stack.clear();
		visited[start] = true;
		stack.push_back(start);

		while (!stack.empty()) {
			const uint idx = stack.back();
			stack.pop_back();
			component.push_back(idx);

			forEachCrossNeighbor(idx, grid, [&](uint neighborIdx) {
				if (!visited[neighborIdx] &&
					!blocked[neighborIdx] &&
					cells[neighborIdx].hasHit) {
					visited[neighborIdx] = true;
					stack.push_back(neighborIdx);
				}
				return true;
			});
		}

		if (component.size() > largestComponent.size()) {
			largestComponent = component;
		}
	}

	if (largestComponent.empty()) {
		return result;
	}

	std::vector<char> localKeepMask(cells.size(), false);
	for (uint idx : largestComponent) {
		localKeepMask[idx] = true;
	}

	for (uint idx : largestComponent) {
		result.area += grid.sideU * grid.sideV;

		const uint row = idx / grid.cols;
		const uint col = idx % grid.cols;

		if (col == 0 || !localKeepMask[idx - 1]) {
			result.perimeter += grid.sideV;
		}
		if (col + 1 == grid.cols || !localKeepMask[idx + 1]) {
			result.perimeter += grid.sideV;
		}
		if (row == 0 || !localKeepMask[idx - grid.cols]) {
			result.perimeter += grid.sideU;
		}
		if (row + 1 == grid.rows || !localKeepMask[idx + grid.cols]) {
			result.perimeter += grid.sideU;
		}
	}

	result.compactness =
		(result.perimeter > 0.0) ?
			(4.0 * M_PI * result.area) /
				(result.perimeter * result.perimeter) :
			0.0;

	if (keepMask != nullptr) {
		*keepMask = std::move(localKeepMask);
	}

	return result;
}

static std::vector<CellData> applyLineCutAndKeepLargest(
	const std::vector<CellData>& cells,
	const GridChoice& grid,
	const std::vector<vcl::uint>& line)
{
	std::vector<char> blocked(cells.size(), false);
	for (vcl::uint idx : line) {
		blocked[idx] = true;
	}

	std::vector<char> keepMask(cells.size(), false);
	largestUnblockedHitComponentShape(cells, grid, blocked, &keepMask);

	std::vector<CellData> result = cells;
	for (vcl::uint idx = 0; idx < result.size(); ++idx) {
		if (!keepMask[idx] && result[idx].hasHit) {
			result[idx].hasHit = false;
			result[idx].isReduced = true;
		}
	}

	return result;
}

struct ChordCutCandidate
{
	vcl::uint a = 0;
	vcl::uint b = 0;
	double length = 0.0;
};

static std::vector<CellData> cutProtrusions(
	std::vector<CellData> cells,
	const GridChoice& grid,
	double maxDistance,
	std::vector<std::vector<vcl::uint>> connectedNeighbors)
{
	using namespace vcl;

	if (cells.size() != grid.rows * grid.cols) {
		return cells;
	}

	for (CellData& cell : cells) {
		cell.hasHit = cell.hasHit && !cell.hasClampedHit;
	}
	
	const double minChordLength =  0.1 * maxDistance;
	const double maxChordLength = 0.25 * maxDistance; ;
	const double maxRemovedAreaRatio = 0.04;
	const double minCompactnessGain = 0.006;
	const double minScore = 0.002;
	const uint maxCuts = 5;
	const uint maxBoundarySamples = 120;
	const uint maxEvaluatedCandidates = 650;

	for (uint cutIteration = 0; cutIteration < maxCuts; ++cutIteration) {
		const HitCellShapeData beforeShape = hitCellShape(cells, grid);
		if (beforeShape.area <= 0.0 || beforeShape.perimeter <= 0.0) {
			break;
		}
		
		//get boundary cells
		const std::vector<uint> boundary = boundaryHitCells(cells, grid);
		if (boundary.size() < 2) {
			break;
		}

		const uint boundaryStride =
			std::max<uint>(
				1,
				static_cast<uint>(
					std::ceil(
						static_cast<double>(boundary.size()) /
						static_cast<double>(maxBoundarySamples))));
		std::vector<uint> sampledBoundary;
		sampledBoundary.reserve(
			(boundary.size() + boundaryStride - 1) / boundaryStride);
		//sample boundary cells
		for (uint i = 0; i < boundary.size(); i += boundaryStride) {
			sampledBoundary.push_back(boundary[i]);
		}

		std::vector<ChordCutCandidate> candidates;
		candidates.reserve(
			sampledBoundary.size() * (sampledBoundary.size() - 1) / 2);
		//for each pair of sampled cells, save possible chords to cut
		//chords are bounded by a min and max length
		for (uint i = 0; i < sampledBoundary.size(); ++i) {
			for (uint j = i + 1; j < sampledBoundary.size(); ++j) {
				const double chordLength =
					gridCellCenterDistance(sampledBoundary[i], sampledBoundary[j], grid);
				if (chordLength < minChordLength ||
					chordLength > maxChordLength) {
					continue;
				}

				candidates.push_back(
					{sampledBoundary[i], sampledBoundary[j], chordLength});
			}
		}

		std::sort(
			candidates.begin(),
			candidates.end(),
			[](const ChordCutCandidate& a, const ChordCutCandidate& b) {
				return a.length < b.length;
			});

		double bestScore = minScore;
		std::vector<uint> bestLine;
		uint evaluatedCandidates = 0;
		std::vector<char> blocked(cells.size(), false);

		for (const ChordCutCandidate& candidate : candidates) {
			if (evaluatedCandidates >= maxEvaluatedCandidates) {
				break;
			}
			//get line between candidate cells
			const std::vector<uint> line =
				rasterizedGridSegment(candidate.a, candidate.b, grid);
			//avoid lines with !cell.hasHit
			if (!lineInsideHitMask(cells, line)) {
				continue;
			}

			++evaluatedCandidates;

			std::fill(blocked.begin(), blocked.end(), false);
			for (uint idx : line) {
				blocked[idx] = true;
			}

			const HitCellShapeData afterShape =
				largestUnblockedHitComponentShape(cells, grid, blocked);
			if (afterShape.area <= 0.0 || afterShape.area >= beforeShape.area) {
				continue;
			}
			
			//reject if it removes too much
			const double removedAreaRatio =
				(beforeShape.area - afterShape.area) / beforeShape.area;
			if (removedAreaRatio <= 0.0 ||
				removedAreaRatio > maxRemovedAreaRatio) {
				continue;
			}

			//reject if it does not improve compactness enough
			const double compactnessGain =
				afterShape.compactness - beforeShape.compactness;
			if (compactnessGain < minCompactnessGain) {
				continue;
			}

			const double perimeterGainRatio =
				(beforeShape.perimeter - afterShape.perimeter) /
				beforeShape.perimeter;
			const double score =
				compactnessGain +
				0.25 * perimeterGainRatio -
				removedAreaRatio;

			if (score > bestScore) {
				bestScore = score;
				bestLine = line;
			}
		}

		if (bestLine.empty()) {
			break;
		}
		
		//apply cut
		cells = applyLineCutAndKeepLargest(cells, grid, bestLine);
	}

	return keepLargestHitComponent(cells, connectedNeighbors);
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
				cells[idx].isReduced = true;
				cells[idx].distance = maxDistance;
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
		if (!candidateCells[idx].hasHit && reducedCells[idx].hasHit) {
			reducedCells[idx].hasHit = false;
			reducedCells[idx].hasClampedHit = false;
			reducedCells[idx].isReduced = true;
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


#endif