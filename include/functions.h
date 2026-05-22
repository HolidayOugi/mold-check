#ifndef VCL_TEST_EXTERNAL_888_MOLD_CHECK_FUNCTIONS_H
#define VCL_TEST_EXTERNAL_888_MOLD_CHECK_FUNCTIONS_H

#include "helper.h"
#include "struct.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
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

	for (uint j = 0; j < cells.size(); ++j) {
		if (i == j || !cells[j].hasHit) {
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

/*

static std::vector<CellData> makeDepthCells(
	const std::vector<CellData>& cells,
	const vcl::Point3d& direction,
	const GridChoice& grid)
{
	if (cells.size() != grid.rows * grid.cols) {
		return {};
	}

	std::vector<CellData> depthCells = cells;
	double distanceSum = 0.0;
	vcl::uint hitCount = 0;

	for (vcl::uint i = 0; i < cells.size(); ++i) {
		if (cells[i].hasHit) {
			depthCells[i].distance = cells[i].clampedDistance;
			depthCells[i].hitPoints = {
				depthCells[i].cellCenter +
				direction * depthCells[i].distance};
			depthCells[i].hasHit = true;

			const vcl::uint row = i / grid.cols;
			const vcl::uint col = i % grid.cols;
			bool hasMissingNeighbor = false;

			if (col > 0 && !cells[i - 1].hasHit) {
				hasMissingNeighbor = true;
			}
			if (col + 1 < grid.cols && !cells[i + 1].hasHit) {
				hasMissingNeighbor = true;
			}
			if (row > 0 && !cells[i - grid.cols].hasHit) {
				hasMissingNeighbor = true;
			}
			if (row + 1 < grid.rows && !cells[i + grid.cols].hasHit) {
				hasMissingNeighbor = true;
			}

			if (hasMissingNeighbor) {
				distanceSum += depthCells[i].distance;
				++hitCount;
			}
		}
		else {
			depthCells[i].hasHit = false;
		}
	}

	if (hitCount == 0) {
		return depthCells;
	}

	const double averageDistance = distanceSum / hitCount;

	for (CellData& depthCell : depthCells) {
		if (!depthCell.hasHit) {
			depthCell.distance = averageDistance;
			depthCell.hitPoints = {
				depthCell.cellCenter + direction * depthCell.distance};
		}
	}

	return depthCells;
}

*/

static std::vector<CellData> reducePoints(
	std::vector<CellData> cells,
	const GridChoice& grid,
	double distanceThreshold = std::numeric_limits<double>::infinity())
{
	using namespace vcl;

	if (cells.size() != grid.rows * grid.cols) {
		return cells;
	}

	std::vector<CellData> candidateCells = cells;
	for (CellData& cell : candidateCells) {
		if (cell.hasClampedHit) {
			cell.hasHit = false;
		}
	}

	erodeHitMaskOnce(candidateCells, grid);
	erodeHitMaskOnce(candidateCells, grid);
	dilateHitMaskOnce(candidateCells, grid);
	dilateHitMaskOnce(candidateCells, grid);

	const std::vector<std::vector<uint>> connectedNeighbors =
		removeDistanceJumpPoints(
			candidateCells,
			grid,
			distanceThreshold);

	candidateCells = keepLargestHitComponent(candidateCells, connectedNeighbors);

	for (uint idx = 0; idx < cells.size(); ++idx) {
		if (!cells[idx].hasClampedHit) {
			cells[idx].hasHit = cells[idx].hasHit && candidateCells[idx].hasHit;
		}
	}

	return cells;
}

/*
static std::vector<CellData> smoothMissingDepthCells(
	const std::vector<CellData>& depthCells,
	const std::vector<CellData>& cells,
	const vcl::Point3d& direction,
	const GridChoice& grid,
	vcl::uint squareSize,
	vcl::uint maxIterations)
{
	using namespace vcl;

	if (depthCells.size() != cells.size() ||
		depthCells.size() != grid.rows * grid.cols) {
		return {};
	}

	std::vector<uint> allCells(depthCells.size());
	std::iota(allCells.begin(), allCells.end(), 0);

	const int radius = static_cast<int>(squareSize / 2);
	std::vector<std::vector<uint>> neighborsByCell(depthCells.size());

	parallelFor(allCells, [&](uint idx) {
		const uint centerRow = idx / grid.cols;
		const uint centerCol = idx % grid.cols;

		for (int rowOffset = -radius; rowOffset <= radius; ++rowOffset) {
			for (int colOffset = -radius; colOffset <= radius; ++colOffset) {
				const int row = static_cast<int>(centerRow) + rowOffset;
				const int col = static_cast<int>(centerCol) + colOffset;

				if (row < 0 ||
					col < 0 ||
					row >= static_cast<int>(grid.rows) ||
					col >= static_cast<int>(grid.cols)) {
					continue;
				}

				const uint neighborIdx =
					static_cast<uint>(row) * grid.cols +
					static_cast<uint>(col);

				if (neighborIdx != idx) {
					neighborsByCell[idx].push_back(neighborIdx);
				}
			}
		}
	});

	std::vector<double> currentDistances(depthCells.size(), 0.0);
	for (uint idx = 0; idx < depthCells.size(); ++idx) {
		currentDistances[idx] = depthCells[idx].distance;
	}

	for (uint iteration = 0; iteration < maxIterations; ++iteration) {
		std::vector<double> nextDistances = currentDistances;

		parallelFor(allCells, [&](uint idx) {
			const bool isMissingDepthCell = !depthCells[idx].hasHit;

			if (!isMissingDepthCell) {
				return;
			}

			double distanceSum = 0.0;
			uint distanceCount = 0;
			double maxAllowedDistance = std::numeric_limits<double>::infinity();

			if (isMissingDepthCell) {
				maxAllowedDistance = cells[idx].clampedDistance;
			}

			for (uint neighborIdx : neighborsByCell[idx]) {
				distanceSum += currentDistances[neighborIdx];
				++distanceCount;

				if (cells[neighborIdx].hasHit &&
					cells[neighborIdx].hasClampedHit) {
					maxAllowedDistance =
						std::min(maxAllowedDistance, cells[neighborIdx].clampedDistance);
				}
			}

			if (distanceCount == 0) {
				return;
			}

			double newDistance = distanceSum / distanceCount;

			if (maxAllowedDistance < std::numeric_limits<double>::infinity()) {
				newDistance = std::min(newDistance, maxAllowedDistance);
			}

			nextDistances[idx] = newDistance;
		});

		currentDistances = std::move(nextDistances);
	}

	std::vector<CellData> result = depthCells;
	for (uint idx = 0; idx < result.size(); ++idx) {
		result[idx].distance = currentDistances[idx];
		result[idx].hitPoints = {
			result[idx].cellCenter + direction * result[idx].distance};
	}

	return result;
}

*/

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

static std::vector<CellData> makeDepthCellsPullPush(
	const std::vector<CellData>& cells,
	const vcl::Point3d& direction,
	const GridChoice& grid,
	double coneCosThreshold,
	float eps)
{
	using namespace vcl;

	if (cells.size() != grid.rows * grid.cols) {
		return {};
	}

	struct PullPushLevel
	{
		uint rows = 0;
		uint cols = 0;
		std::vector<double> distances;
		std::vector<double> weights;
	};

	std::vector<CellData> depthCells = cells;
	std::vector<PullPushLevel> pyramid;
	pyramid.push_back({
		grid.rows,
		grid.cols,
		std::vector<double>(cells.size(), 0.0),
		std::vector<double>(cells.size(), 0.0)});

	uint hitCount = 0;
	for (uint idx = 0; idx < depthCells.size(); ++idx) {
		depthCells[idx].distance = cells[idx].clampedDistance;
		depthCells[idx].hitPoints = {
			depthCells[idx].cellCenter + direction * depthCells[idx].distance};
		if (cells[idx].hasHit) {
			pyramid[0].distances[idx] = depthCells[idx].distance;
			pyramid[0].weights[idx] = 1.0;
			++hitCount;
		}
	}

	if (hitCount == 0) {
		return cells;
	}

	while (pyramid.back().rows > 1 || pyramid.back().cols > 1) {
		const PullPushLevel& fine = pyramid.back();
		const uint coarseRows = (fine.rows + 1) / 2;
		const uint coarseCols = (fine.cols + 1) / 2;

		PullPushLevel coarse{
			coarseRows,
			coarseCols,
			std::vector<double>(coarseRows * coarseCols, 0.0),
			std::vector<double>(coarseRows * coarseCols, 0.0)};

		for (uint row = 0; row < coarseRows; ++row) {
			for (uint col = 0; col < coarseCols; ++col) {
				double weightedDistanceSum = 0.0;
				double weightSum = 0.0;

				for (uint rowOffset = 0; rowOffset < 2; ++rowOffset) {
					for (uint colOffset = 0; colOffset < 2; ++colOffset) {
						const uint fineRow = row * 2 + rowOffset;
						const uint fineCol = col * 2 + colOffset;

						if (fineRow >= fine.rows || fineCol >= fine.cols) {
							continue;
						}

						const uint fineIdx = fineRow * fine.cols + fineCol;
						const double weight = fine.weights[fineIdx];

						if (weight <= 0.0) {
							continue;
						}

						weightedDistanceSum += fine.distances[fineIdx] * weight;
						weightSum += weight;
					}
				}

				const uint coarseIdx = row * coarseCols + col;
				if (weightSum > 0.0) {
					coarse.distances[coarseIdx] = weightedDistanceSum / weightSum;
					coarse.weights[coarseIdx] = weightSum;
				}
			}
		}

		pyramid.push_back(std::move(coarse));
	}

	for (uint levelIdx = static_cast<uint>(pyramid.size() - 1);
		 levelIdx > 0;
		 --levelIdx) {
		const PullPushLevel& coarse = pyramid[levelIdx];
		PullPushLevel& fine = pyramid[levelIdx - 1];

		for (uint row = 0; row < fine.rows; ++row) {
			for (uint col = 0; col < fine.cols; ++col) {
				const uint idx = row * fine.cols + col;

				if (levelIdx == 1) {
					if (cells[idx].hasHit && !cells[idx].hasClampedHit) {
						fine.distances[idx] = cells[idx].clampedDistance;
						fine.weights[idx] = 1.0;
						continue;
					}
				}
				else if (fine.weights[idx] > 0.0) {
					continue;
				}

				const double coarseRow =
					(static_cast<double>(row) + 0.5) *
					static_cast<double>(coarse.rows) /
					static_cast<double>(fine.rows) - 0.5;
				const double coarseCol =
					(static_cast<double>(col) + 0.5) *
					static_cast<double>(coarse.cols) /
					static_cast<double>(fine.cols) - 0.5;
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

						const uint sampleIdx =
							static_cast<uint>(sampleRow) * coarse.cols +
							static_cast<uint>(sampleCol);
						const double interpolationWeight =
							std::max(
								0.0,
								1.0 - std::abs(coarseRow - sampleRow)) *
							std::max(
								0.0,
								1.0 - std::abs(coarseCol - sampleCol));
						const double weight =
							interpolationWeight * coarse.weights[sampleIdx];

						weightedDistanceSum += coarse.distances[sampleIdx] * weight;
						weightSum += weight;
					}
				}

				if (weightSum > 0.0) {
					fine.distances[idx] = weightedDistanceSum / weightSum;
				}
				else {
					const uint parentRow = std::min(row / 2, coarse.rows - 1);
					const uint parentCol = std::min(col / 2, coarse.cols - 1);
					fine.distances[idx] =
						coarse.distances[parentRow * coarse.cols + parentCol];
				}

				if (levelIdx == 1 && cells[idx].hasClampedHit) {
					fine.distances[idx] =
						std::min(fine.distances[idx], cells[idx].clampedDistance);
				}
				fine.weights[idx] = 1.0;
			}
		}
	}

	for (uint idx = 0; idx < depthCells.size(); ++idx) {
		depthCells[idx].distance = pyramid[0].distances[idx];
		if (cells[idx].hasClampedHit) {
			depthCells[idx].distance =
				std::min(depthCells[idx].distance, cells[idx].clampedDistance);
		}
		depthCells[idx].hitPoints = {
			depthCells[idx].cellCenter + direction * depthCells[idx].distance};
		depthCells[idx].hasHit = cells[idx].hasHit;
	}

	depthCells =
		fixDepthCellConeViolations(
			depthCells,
			direction,
			coneCosThreshold,
			eps);

	for (uint idx = 0; idx < depthCells.size(); ++idx) {
		if (depthCells[idx].hasClampedHit) {
			depthCells[idx].distance =
				std::min(depthCells[idx].distance, depthCells[idx].clampedDistance);
		}
		depthCells[idx].hitPoints = {
			depthCells[idx].cellCenter + direction * depthCells[idx].distance};
		depthCells[idx].hasHit = cells[idx].hasHit;
	}

	return depthCells;
}

#endif
