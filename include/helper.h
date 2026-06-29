#ifndef VCL_TEST_EXTERNAL_888_MOLD_CHECK_HELPER_H
#define VCL_TEST_EXTERNAL_888_MOLD_CHECK_HELPER_H

#include "struct.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <unordered_set>
#include <vector>

#include <vclib/base/parallel.h>
#include <vclib/meshes.h>

static vcl::Point3d computeHitPoint(
	const vcl::PolyMesh& m,
	vcl::uint faceId,
	vcl::uint triId,
	const vcl::Point3f& baryCoords,
	const vcl::Point3d& invalidPoint)
{
	using namespace vcl;

	if (faceId == UINT_NULL) {
		return invalidPoint;
	}

	const auto& face = m.face(faceId);
	const std::vector<uint> faceTriangulation = earCut(face);

	if (triId * 3 + 2 >= faceTriangulation.size()) {
		return invalidPoint;
	}

	const uint vi0 = faceTriangulation[triId * 3 + 0];
	const uint vi1 = faceTriangulation[triId * 3 + 1];
	const uint vi2 = faceTriangulation[triId * 3 + 2];

	const Point3d& p0 = face.vertex(vi0)->position();
	const Point3d& p1 = face.vertex(vi1)->position();
	const Point3d& p2 = face.vertex(vi2)->position();

	return
		p0 * baryCoords.x() +
		p1 * baryCoords.y() +
		p2 * baryCoords.z();
}

static bool isWithinPlaneAngle(
	const vcl::Point3d& point,
	const vcl::Point3d& other,
	const vcl::Point3d& direction,
	double coneCosThreshold,
	float eps)
{	

	using namespace vcl;

	const Point3d directionToOther = other - point;
	const double directionToOtherNorm = directionToOther.norm();

	if (directionToOtherNorm < eps) {
		return true;
	}

	const Point3d dirToPlane = -direction;

	const Point3d dirToOtherNormalized =
		directionToOther / directionToOtherNorm;

	const double cosBetween =
		dirToOtherNormalized.dot(dirToPlane);

	return cosBetween > coneCosThreshold - eps;
}

static double coneBoundaryStep(
	const vcl::Point3d& a,
	const vcl::Point3d& b,
	const vcl::Point3d& direction,
	double coneCosThreshold,
	float eps)
{

	using namespace vcl;

	const Point3d ab = b - a;
	const double abNorm = ab.norm();

	if (abNorm < eps) {
		return 0.0;
	}

	// If b is already outside the cone of a, no displacement is needed.
	// The cone has its apex at a and opens toward the plane, i.e. along -direction.
	const double cosAngle = (ab / abNorm).dot(-direction);

	if (cosAngle <= coneCosThreshold) {
		return 0.0;
	}

	// We want to move a toward the plane:
	//
	//     newA = a - t * direction
	//
	// such that b lies exactly on the cone boundary of newA.
	//
	// The cone boundary condition is:
	//
	//     dot((b - newA) / |b - newA|, -direction) = coneCosThreshold
	//
	// Since:
	//
	//     b - newA = b - (a - t * direction)
	//              = (b - a) + t * direction
	//
	// let:
	//
	//     v = b - a
	//     d = direction
	//
	// so the displaced vector becomes:
	//
	//     v(t) = v + t * d
	//
	// We solve:
	//
	//     dot(v(t) / |v(t)|, -d) = coneCosThreshold
	//
	// which expands to:
	//
	//     -(dot(v, d) + t) / sqrt(|v|^2 + 2t dot(v, d) + t^2)
	//         = coneCosThreshold
	//
	// Squaring both sides gives a quadratic equation in t.

	const Point3d vec = b - a;

	const double vecDotDirection = vec.dot(direction);
	const double vecNorm2 = vec.dot(vec);

	const double cos2 = coneCosThreshold * coneCosThreshold;
	const double sin2 = 1.0 - cos2;

	// Quadratic equation:
	//
	//     sin^2(theta) * t^2
	//   + 2 * dot(v, d) * sin^2(theta) * t
	//   + dot(v, d)^2 - cos^2(theta) * |v|^2
	//   = 0
	//
	// where cos(theta) = coneCosThreshold.

	const double aEq = sin2;
	const double bEq = 2.0 * vecDotDirection * sin2;
	const double cEq =
		vecDotDirection * vecDotDirection -
		cos2 * vecNorm2;

	const double discriminant =
		bEq * bEq - 4.0 * aEq * cEq;

	if (discriminant <= 0.0) {
		// No valid real solution was found.
		// Fall back to the maximum safe displacement along the direction axis.
		return std::abs(vecDotDirection);
	}

	const double sqrtDisc = std::sqrt(discriminant);

	const double t1 =
		(-bEq - sqrtDisc) / (2.0 * aEq);

	const double t2 =
		(-bEq + sqrtDisc) / (2.0 * aEq);

	// Pick the smallest positive displacement.
	double t = std::numeric_limits<double>::max();

	if (t1 > eps) {
		t = std::min(t, t1);
	}

	if (t2 > eps) {
		t = std::min(t, t2);
	}

	if (t == std::numeric_limits<double>::max()) {
		return 0.0;
	}

	// Do not move farther than the projection of b - a along direction.
	// This prevents overshooting past the plane-aligned limit.
	if (t > std::abs(vecDotDirection)) {
		return std::abs(vecDotDirection);
	}

	return t;
}

static HitCellShapeData hitCellShape(
	const std::vector<CellData>& cells,
	const GridChoice& grid)
{
	using namespace vcl;

	HitCellShapeData result;

	if (cells.size() != grid.rows * grid.cols) {
		return result;
	}

	for (uint idx = 0; idx < cells.size(); ++idx) {
		if (!cells[idx].hasHit) {
			continue;
		}

		result.area += grid.sideU * grid.sideV;

		const uint row = idx / grid.cols;
		const uint col = idx % grid.cols;

		if (col == 0 || !cells[idx - 1].hasHit) {
			result.perimeter += grid.sideV;
		}
		if (col + 1 == grid.cols || !cells[idx + 1].hasHit) {
			result.perimeter += grid.sideV;
		}
		if (row == 0 || !cells[idx - grid.cols].hasHit) {
			result.perimeter += grid.sideU;
		}
		if (row + 1 == grid.rows || !cells[idx + grid.cols].hasHit) {
			result.perimeter += grid.sideU;
		}
	}

	result.compactness =
		(result.perimeter > 0.0) ?
			(4.0 * M_PI * result.area) / (result.perimeter * result.perimeter) :
			0.0;

	return result;
}

static std::vector<vcl::uint> crossNeighborIndices(
	vcl::uint idx,
	const GridChoice& grid)
{
	using namespace vcl;

	const uint row = idx / grid.cols;
	const uint col = idx % grid.cols;

	std::vector<uint> neighbors;
	neighbors.reserve(4);

	if (col > 0) {
		neighbors.push_back(idx - 1);
	}
	if (col + 1 < grid.cols) {
		neighbors.push_back(idx + 1);
	}
	if (row > 0) {
		neighbors.push_back(idx - grid.cols);
	}
	if (row + 1 < grid.rows) {
		neighbors.push_back(idx + grid.cols);
	}

	return neighbors;
}

static std::vector<vcl::uint> squareNeighborIndices(
	vcl::uint idx,
	const GridChoice& grid,
	vcl::uint squareSize)
{
	using namespace vcl;

	const uint centerRow = idx / grid.cols;
	const uint centerCol = idx % grid.cols;
	const int radius = static_cast<int>(squareSize / 2);

	std::vector<uint> neighbors;
	neighbors.reserve(squareSize * squareSize - 1);

	for (int rowOffset = -radius; rowOffset <= radius; ++rowOffset) {
		for (int colOffset = -radius; colOffset <= radius; ++colOffset) {
			if (rowOffset == 0 && colOffset == 0) {
				continue;
			}

			const int row = static_cast<int>(centerRow) + rowOffset;
			const int col = static_cast<int>(centerCol) + colOffset;

			if (row < 0 ||
				col < 0 ||
				row >= static_cast<int>(grid.rows) ||
				col >= static_cast<int>(grid.cols)) {
				continue;
			}

			neighbors.push_back(
				static_cast<uint>(row) * grid.cols +
				static_cast<uint>(col));
		}
	}

	return neighbors;
}

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
		for (uint neighborIdx : crossNeighborIndices(idx, grid)) {
			if (!cells[neighborIdx].hasHit) {
				keep = false;
				break;
			}
		}

		nextHasHit[idx] = keep;
	});

	parallelFor(allCells, [&](uint idx) {
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
		for (uint neighborIdx : crossNeighborIndices(idx, grid)) {
			if (cells[neighborIdx].hasHit) {
				add = true;
				break;
			}
		}

		nextHasHit[idx] = add;
	});

	parallelFor(allCells, [&](uint idx) {
		cells[idx].hasHit = nextHasHit[idx];
	});
}

static std::vector<std::vector<vcl::uint>> removeDistanceJumpPoints(
	const std::vector<CellData>& cells,
	const GridChoice& grid,
	double distanceThreshold)
{
	using namespace vcl;

	if (distanceThreshold <= 0.0) {
		distanceThreshold = std::numeric_limits<double>::infinity();
	}

	std::vector<std::vector<uint>> connectedNeighbors(cells.size());
	for (uint idx = 0; idx < cells.size(); ++idx) {
		if (!cells[idx].hasHit) {
			continue;
		}

		for (uint neighborIdx : crossNeighborIndices(idx, grid)) {
			if (neighborIdx <= idx || !cells[neighborIdx].hasHit) {
				continue;
			}

			if (std::abs(cells[idx].distance - cells[neighborIdx].distance) <=
				distanceThreshold) {
				connectedNeighbors[idx].push_back(neighborIdx);
				connectedNeighbors[neighborIdx].push_back(idx);
			}
		}
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

	std::unordered_set<uint> keep(
		largestComponent.begin(), largestComponent.end());
	std::vector<CellData> result = cells;

	for (uint idx = 0; idx < result.size(); ++idx) {
		if (result[idx].hasHit && keep.count(idx) == 0) {
			result[idx].hasHit = false;
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

	for (uint neighborIdx : crossNeighborIndices(idx, grid)) {
		if (!cells[neighborIdx].hasHit) {
			return true;
		}
	}

	return false;
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

			for (uint neighborIdx : crossNeighborIndices(idx, grid)) {
				if (!visited[neighborIdx] &&
					!blocked[neighborIdx] &&
					cells[neighborIdx].hasHit) {
					visited[neighborIdx] = true;
					stack.push_back(neighborIdx);
				}
			}
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
		if (!keepMask[idx]) {
			result[idx].hasHit = false;
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
		//sample boundary cells
		for (uint i = 0; i < boundary.size(); i += boundaryStride) {
			sampledBoundary.push_back(boundary[i]);
		}

		std::vector<ChordCutCandidate> candidates;
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

static double interpolateFromCoarseLevel(
	const PullPushLevel& coarse,
	vcl::uint fineRow,
	vcl::uint fineCol,
	vcl::uint fineRows,
	vcl::uint fineCols)
{
	const double coarseRow =
		(static_cast<double>(fineRow) + 0.5) *
		static_cast<double>(coarse.rows) /
		static_cast<double>(fineRows) - 0.5;

	const double coarseCol =
		(static_cast<double>(fineCol) + 0.5) *
		static_cast<double>(coarse.cols) /
		static_cast<double>(fineCols) - 0.5;

	const int row0 = static_cast<int>(std::floor(coarseRow));
	const int col0 = static_cast<int>(std::floor(coarseCol));

	double weightedDistanceSum = 0.0;
	double weightSum = 0.0;

	for (int rowOffset = 0; rowOffset <= 1; ++rowOffset) {
		for (int colOffset = 0; colOffset <= 1; ++colOffset) {
			const int sampleRow = row0 + rowOffset;
			const int sampleCol = col0 + colOffset;

			if (sampleRow < 0 ||
				sampleCol < 0 ||
				sampleRow >= static_cast<int>(coarse.rows) ||
				sampleCol >= static_cast<int>(coarse.cols)) {
				continue;
			}

			const double rowWeight =
				std::max(0.0, 1.0 - std::abs(coarseRow - sampleRow));

			const double colWeight =
				std::max(0.0, 1.0 - std::abs(coarseCol - sampleCol));

			const double interpolationWeight = rowWeight * colWeight;

			const vcl::uint sampleIndex =
				static_cast<vcl::uint>(sampleRow) * coarse.cols +
				static_cast<vcl::uint>(sampleCol);

			const double weight =
				interpolationWeight * coarse.weights[sampleIndex];

			weightedDistanceSum += coarse.distances[sampleIndex] * weight;
			weightSum += weight;
		}
	}

	if (weightSum > 0.0) {
		return weightedDistanceSum / weightSum;
	}

	const vcl::uint parentRow = std::min(fineRow / 2, coarse.rows - 1);
	const vcl::uint parentCol = std::min(fineCol / 2, coarse.cols - 1);

	return coarse.distances[parentRow * coarse.cols + parentCol];
}

#endif
