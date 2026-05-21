#ifndef VCL_TEST_EXTERNAL_888_MOLD_CHECK_HELPER_H
#define VCL_TEST_EXTERNAL_888_MOLD_CHECK_HELPER_H

#include "struct.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <vector>

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

static bool isSameDistanceCell(
	const std::vector<CellData>& cells,
	const std::vector<CellData>& clampedCells,
	vcl::uint idx,
	float eps)
{
	return cells[idx].hasHit &&
		   clampedCells[idx].hasHit &&
		   std::abs(cells[idx].distance - clampedCells[idx].distance) <= eps;
}

static void pushNeighbor(
	std::vector<vcl::uint>& stack,
	std::vector<bool>& visited,
	const std::vector<CellData>& cells,
	const std::vector<CellData>& clampedCells,
	vcl::uint neighbor,
	float eps)
{
	if (!visited[neighbor] &&
		isSameDistanceCell(cells, clampedCells, neighbor, eps)) {
		visited[neighbor] = true;
		stack.push_back(neighbor);
	}
}

static double componentGridPerimeter(
	const std::vector<vcl::uint>& componentIndices,
	const GridChoice& grid)
{
	using namespace vcl;

	double perimeter = 0.0;
	std::unordered_set<uint> componentSet(
		componentIndices.begin(), componentIndices.end());

	for (uint idx : componentIndices) {
		const uint row = idx / grid.cols;
		const uint col = idx % grid.cols;

		if (col == 0 || componentSet.count(idx - 1) == 0) {
			perimeter += grid.sideV;
		}
		if (col + 1 == grid.cols || componentSet.count(idx + 1) == 0) {
			perimeter += grid.sideV;
		}
		if (row == 0 || componentSet.count(idx - grid.cols) == 0) {
			perimeter += grid.sideU;
		}
		if (row + 1 == grid.rows || componentSet.count(idx + grid.cols) == 0) {
			perimeter += grid.sideU;
		}
	}

	return perimeter;
}

static std::vector<vcl::uint> squareNeighborIndices(
	vcl::uint idx,
	const GridChoice& grid,
	vcl::uint squareSize)
{
	using namespace vcl;

	const uint row = idx / grid.cols;
	const uint col = idx % grid.cols;
	const int radius = static_cast<int>(squareSize / 2);

	std::vector<uint> neighbors;
	neighbors.reserve(squareSize * squareSize);

	for (int rowOffset = -radius; rowOffset <= radius; ++rowOffset) {
		for (int colOffset = -radius; colOffset <= radius; ++colOffset) {
			const int neighborRow = static_cast<int>(row) + rowOffset;
			const int neighborCol = static_cast<int>(col) + colOffset;

			if (neighborRow < 0 ||
				neighborCol < 0 ||
				neighborRow >= static_cast<int>(grid.rows) ||
				neighborCol >= static_cast<int>(grid.cols)) {
				continue;
			}

			const uint neighborIdx =
				static_cast<uint>(neighborRow) * grid.cols +
				static_cast<uint>(neighborCol);

			if (neighborIdx == idx) {
				continue;
			}

			neighbors.push_back(neighborIdx);
		}
	}

	return neighbors;
}

static void erodeHitMaskOnce(
	std::vector<CellData>& cells,
	const GridChoice& grid,
	vcl::uint squareSize)
{
	using namespace vcl;

	std::vector<bool> nextHasHit(cells.size(), false);

	for (uint idx = 0; idx < cells.size(); ++idx) {
		if (!cells[idx].hasHit) {
			continue;
		}

		bool keep = true;
		for (uint neighborIdx : squareNeighborIndices(idx, grid, squareSize)) {
			if (!cells[neighborIdx].hasHit) {
				keep = false;
				break;
			}
		}

		nextHasHit[idx] = keep;
	}

	for (uint idx = 0; idx < cells.size(); ++idx) {
		cells[idx].hasHit = nextHasHit[idx];
	}
}

static void dilateHitMaskOnce(
	std::vector<CellData>& cells,
	const GridChoice& grid,
	vcl::uint squareSize)
{
	using namespace vcl;

	std::vector<bool> nextHasHit(cells.size(), false);

	for (uint idx = 0; idx < cells.size(); ++idx) {
		bool add = cells[idx].hasHit;
		for (uint neighborIdx : squareNeighborIndices(idx, grid, squareSize)) {
			if (cells[neighborIdx].hasHit) {
				add = true;
				break;
			}
		}

		nextHasHit[idx] = add;
	}

	for (uint idx = 0; idx < cells.size(); ++idx) {
		cells[idx].hasHit = nextHasHit[idx];
	}
}

static std::vector<CellData> removeDistanceJumpPoints(
	const std::vector<CellData>& cells,
	const GridChoice& grid,
	vcl::uint squareSize,
	double distanceThreshold)
{
	using namespace vcl;

	if (distanceThreshold <= 0.0) {
		return cells;
	}

	std::vector<CellData> result = cells;

	for (uint idx = 0; idx < cells.size(); ++idx) {
		if (!cells[idx].hasHit) {
			continue;
		}

		for (uint neighborIdx : squareNeighborIndices(idx, grid, squareSize)) {
			if (!cells[neighborIdx].hasHit) {
				continue;
			}

			if (std::abs(cells[idx].distance - cells[neighborIdx].distance) >
				distanceThreshold) {
				result[idx].hasHit = false;
				break;
			}
		}
	}

	return result;
}

static std::vector<CellData> keepLargestHitComponent(
	const std::vector<CellData>& cells,
	const GridChoice& grid,
	vcl::uint squareSize)
{
	using namespace vcl;

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

			for (uint neighborIdx : squareNeighborIndices(idx, grid, squareSize)) {
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

#endif
