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

			return result;
		}
	
	}

	CellData result = cell;
	result.distance = maxDistance;
	result.hitPoints = {invalidPoint};
	result.hasHit = false;

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

	return cell;
}

static CellData computeClampedCell(
	vcl::uint i,
	const std::vector<CellData>& cells,
	const vcl::Point3d& planePoint,
	const vcl::Point3d& direction,
	double coneCosThreshold,
	float eps)
{
	using namespace vcl;

	const CellData baseCell = cells[i];

	if (!baseCell.hasHit) {
		return baseCell;
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
		return baseCell;
	}

	const Point3d currentPoint =
		original - direction * requiredT;

	const double distanceToPlane =
		std::abs((currentPoint - planePoint).dot(direction));

	CellData result = baseCell;
	result.distance = distanceToPlane;
	result.hitPoints[0] = currentPoint;
	result.hasHit = true;
	return result;
}

static std::vector<CellData> makeDepthCells(
	const std::vector<CellData>& cells,
	const std::vector<CellData>& clampedCells,
	const vcl::Point3d& direction,
	const GridChoice& grid)
{
	if (cells.size() != clampedCells.size() ||
		cells.size() != grid.rows * grid.cols) {
		return {};
	}

	std::vector<CellData> depthCells = cells;
	double distanceSum = 0.0;
	vcl::uint hitCount = 0;

	for (vcl::uint i = 0; i < cells.size(); ++i) {
		if (cells[i].hasHit) {
			depthCells[i] = clampedCells[i];
			depthCells[i].hitPoints = {
				depthCells[i].cellCenter +
				direction * depthCells[i].distance};

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

static std::vector<CellData> reducePoints(
	std::vector<CellData> cells,
	const GridChoice& grid,
	double distanceThreshold = std::numeric_limits<double>::infinity())
{
	using namespace vcl;

	if (cells.size() != grid.rows * grid.cols) {
		return cells;
	}

	erodeHitMaskOnce(cells, grid);
	erodeHitMaskOnce(cells, grid);
	dilateHitMaskOnce(cells, grid);
	dilateHitMaskOnce(cells, grid);

	const std::vector<std::vector<uint>> connectedNeighbors =
		removeDistanceJumpPoints(
		cells,
		grid,
		distanceThreshold);

	cells = keepLargestHitComponent(cells, connectedNeighbors);

	return cells;
}

static std::vector<CellData> smoothMissingDepthCells(
	const std::vector<CellData>& depthCells,
	const std::vector<CellData>& cells,
	const std::vector<CellData>& clampedCells,
	const vcl::Point3d& direction,
	const GridChoice& grid,
	vcl::uint squareSize,
	vcl::uint maxIterations)
{
	using namespace vcl;

	if (depthCells.size() != cells.size() ||
		depthCells.size() != clampedCells.size() ||
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
			const bool isClampedCloserCell =
				cells[idx].hasHit &&
				clampedCells[idx].hasHit &&
				clampedCells[idx].distance < cells[idx].distance;

			if (!isMissingDepthCell && !isClampedCloserCell) {
				return;
			}

			double distanceSum = 0.0;
			uint distanceCount = 0;
			double maxAllowedDistance = std::numeric_limits<double>::infinity();

			if (isMissingDepthCell) {
				maxAllowedDistance = cells[idx].distance;
			}

			for (uint neighborIdx : neighborsByCell[idx]) {
				distanceSum += currentDistances[neighborIdx];
				++distanceCount;

				if (cells[neighborIdx].hasHit &&
					clampedCells[neighborIdx].hasHit &&
					clampedCells[neighborIdx].distance != cells[neighborIdx].distance) {
					maxAllowedDistance =
						std::min(maxAllowedDistance, clampedCells[neighborIdx].distance);
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
	});

	return depthCells;
}

#endif
