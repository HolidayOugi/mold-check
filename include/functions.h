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

#include <vclib/algorithms/mesh/stat/geometry.h>
#include <vclib/algorithms/mesh/update/bounding_box.h>
#include <vclib/embree/scene.h>
#include <vclib/meshes.h>

static std::tuple<vcl::Point3d, vcl::Point3d> expandedBoundingBox(
	const vcl::PolyMesh& m,
	double marginFactor)
{
	using namespace vcl;

	if (m.vertexCount() == 0) {
		return {Point3d(), Point3d()};
	}

	Point3d min = m.vertices().begin()->position();
	Point3d max = min;

	for (const auto& v : m.vertices()) {
		const Point3d& p = v.position();
		for (uint i = 0; i < 3; ++i) {
			min(i) = std::min(min(i), p(i));
			max(i) = std::max(max(i), p(i));
		}
	}

	const double margin = marginFactor * (max - min).norm();
	min -= Point3d(margin, margin, margin);
	max += Point3d(margin, margin, margin);

	return {min, max};
}

static void addBoxTriangles(
	vcl::PolyMesh& mold,
	const vcl::Point3d& min,
	const vcl::Point3d& max)
{
	using namespace vcl;

	const uint v0 = mold.addVertex(Point3d(min.x(), min.y(), min.z()));
	const uint v1 = mold.addVertex(Point3d(max.x(), min.y(), min.z()));
	const uint v2 = mold.addVertex(Point3d(min.x(), max.y(), min.z()));
	const uint v3 = mold.addVertex(Point3d(max.x(), max.y(), min.z()));
	const uint v4 = mold.addVertex(Point3d(min.x(), min.y(), max.z()));
	const uint v5 = mold.addVertex(Point3d(max.x(), min.y(), max.z()));
	const uint v6 = mold.addVertex(Point3d(min.x(), max.y(), max.z()));
	const uint v7 = mold.addVertex(Point3d(max.x(), max.y(), max.z()));

	mold.addFace(v0, v2, v1);
	mold.addFace(v3, v1, v2);
	mold.addFace(v0, v4, v2);
	mold.addFace(v6, v2, v4);
	mold.addFace(v0, v1, v4);
	mold.addFace(v5, v4, v1);
	mold.addFace(v7, v6, v5);
	mold.addFace(v4, v5, v6);
	mold.addFace(v7, v3, v6);
	mold.addFace(v2, v6, v3);
	mold.addFace(v7, v5, v3);
	mold.addFace(v1, v3, v5);
}

static vcl::PolyMesh squareMold(const vcl::PolyMesh& m, double marginFactor)
{
	using namespace vcl;

	if (m.vertexCount() == 0) {
		return {};
	}

	const auto [min, max] = expandedBoundingBox(m, marginFactor);

	PolyMesh mold;
	addBoxTriangles(mold, min, max);

	updateBoundingBox(mold);
	return mold;
}

static double squareMoldVolume(const vcl::PolyMesh& m, double marginFactor)
{
	using namespace vcl;

	if (m.vertexCount() == 0) {
		return 0.0;
	}

	const auto [min, max] = expandedBoundingBox(m, marginFactor);
	const Point3d boxSize = max - min;
	const double boxVolume = boxSize.x() * boxSize.y() * boxSize.z();

	return boxVolume - std::abs(volume(m));
}

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
		rayHits.erase(
			std::remove_if(
				rayHits.begin(),
				rayHits.end(),
				[&](const auto& rayHit) {
					const uint hitFaceId = std::get<0>(rayHit);
					return !seenFaceIds.insert(hitFaceId).second;
				}),
			rayHits.end());

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
	const std::vector<CellData>& clampedCells)
{
	if (cells.size() != clampedCells.size()) {
		return {};
	}

	std::vector<CellData> depthCells = cells;
	double distanceSum = 0.0;
	vcl::uint hitCount = 0;

	for (vcl::uint i = 0; i < cells.size(); ++i) {
		if (cells[i].hasHit) {
			depthCells[i] = clampedCells[i];
			distanceSum += depthCells[i].distance;
			++hitCount;
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
		}
	}

	return depthCells;
}

static std::vector<CellData> smoothMissingDepthCells(
	const std::vector<CellData>& depthCells,
	const std::vector<CellData>& cells,
	const std::vector<CellData>& clampedCells,
	const GridChoice& grid,
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

	std::vector<CellData> currentDepthCells = depthCells;

	for (uint iteration = 0; iteration < maxIterations; ++iteration) {
		std::vector<CellData> nextDepthCells = currentDepthCells;

		parallelFor(allCells, [&](uint idx) {
			if (currentDepthCells[idx].hasHit) {
				return;
			}

			const uint centerRow = idx / grid.cols;
			const uint centerCol = idx % grid.cols;

			double distanceSum = 0.0;
			uint distanceCount = 0;
			double maxAllowedDistance = std::numeric_limits<double>::infinity();

			const uint minRow = (centerRow > 0) ? centerRow - 1 : centerRow;
			const uint minCol = (centerCol > 0) ? centerCol - 1 : centerCol;
			const uint maxRow = std::min(centerRow + 1, grid.rows - 1);
			const uint maxCol = std::min(centerCol + 1, grid.cols - 1);

			for (uint row = minRow; row <= maxRow; ++row) {
				for (uint col = minCol; col <= maxCol; ++col) {
					const uint neighborIdx = row * grid.cols + col;

					if (neighborIdx == idx) {
						continue;
					}

					distanceSum += currentDepthCells[neighborIdx].distance;
					++distanceCount;

					if (cells[neighborIdx].hasHit &&
						clampedCells[neighborIdx].hasHit &&
						clampedCells[neighborIdx].distance != cells[neighborIdx].distance) {
						maxAllowedDistance =
							std::min(maxAllowedDistance, clampedCells[neighborIdx].distance);
					}
				}
			}

			if (distanceCount == 0) {
				return;
			}

			double newDistance = distanceSum / distanceCount;

			if (maxAllowedDistance < std::numeric_limits<double>::infinity()) {
				newDistance = std::min(newDistance, maxAllowedDistance);
			}

			nextDepthCells[idx].distance = newDistance;
		});

		currentDepthCells = std::move(nextDepthCells);
	}

	return currentDepthCells;
}

static ConnectedComponentData largestConnectedComponent(
	const std::vector<CellData>& cells,
	const std::vector<CellData>& clampedCells,
	const GridChoice& grid,
	float eps)
{
	using namespace vcl;

	ConnectedComponentData result;

	if (cells.size() != clampedCells.size() ||
		cells.size() != grid.rows * grid.cols) {
		return result;
	}

	std::vector<bool> visited(cells.size(), false);

	//BFS

	for (uint start = 0; start < cells.size(); ++start) {
		if (visited[start] ||
			!isSameDistanceCell(cells, clampedCells, start, eps)) {
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

			const uint row = idx / grid.cols;
			const uint col = idx % grid.cols;

			if (col > 0) {
				const uint neighbor = idx - 1;
				pushNeighbor(
					stack, visited, cells, clampedCells, neighbor, eps);
			}
			if (col + 1 < grid.cols) {
				const uint neighbor = idx + 1;
				pushNeighbor(
					stack, visited, cells, clampedCells, neighbor, eps);
			}
			if (row > 0) {
				const uint neighbor = idx - grid.cols;
				pushNeighbor(
					stack, visited, cells, clampedCells, neighbor, eps);
			}
			if (row + 1 < grid.rows) {
				const uint neighbor = idx + grid.cols;
				pushNeighbor(
					stack, visited, cells, clampedCells, neighbor, eps);
			}
		}

		if (component.size() > result.indices.size()) {
			result.indices = std::move(component);
		}
	}

	for (uint idx : result.indices) {
		result.area += grid.sideU * grid.sideV;
	}

	result.perimeter = componentGridPerimeter(result.indices, grid);
	result.compactness =
		(result.perimeter > 0.0) ? (result.area / result.perimeter) : 0.0;

	return result;
}

#endif
