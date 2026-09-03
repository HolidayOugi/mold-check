#ifndef VCL_TEST_EXTERNAL_888_MOLD_CHECK_HELPER_H
#define VCL_TEST_EXTERNAL_888_MOLD_CHECK_HELPER_H

#include "struct.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <numeric>
#include <string>
#include <unordered_set>
#include <vector>

#include <vclib/base/parallel.h>
#include <vclib/io.h>
#include <vclib/meshes.h>

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

template<typename Func>
static bool forEachCrossNeighbor(
	vcl::uint idx,
	const GridChoice& grid,
	Func&& func)
{
	using namespace vcl;

	const uint row = idx / grid.cols;
	const uint col = idx % grid.cols;

	if (col > 0 && !func(idx - 1)) {
		return false;
	}
	if (col + 1 < grid.cols && !func(idx + 1)) {
		return false;
	}
	if (row > 0 && !func(idx - grid.cols)) {
		return false;
	}
	if (row + 1 < grid.rows && !func(idx + grid.cols)) {
		return false;
	}

	return true;
}

template<typename Func>
static bool forEachSquareNeighbor(
	vcl::uint idx,
	const GridChoice& grid,
	vcl::uint squareSize,
	Func&& func)
{
	using namespace vcl;

	const uint centerRow = idx / grid.cols;
	const uint centerCol = idx % grid.cols;
	const int radius = static_cast<int>(squareSize / 2);

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

			if (!func(
					static_cast<uint>(row) * grid.cols +
					static_cast<uint>(col))) {
				return false;
			}
		}
	}

	return true;
}

// Return true when a four-connected region of white cells is fully enclosed
// inside the grid instead of being connected to its outer boundary.
static bool hasEnclosedWhiteHole(
	const std::vector<CellData>& cells,
	const GridChoice& grid,
	const vcl::Point3d& direction,
	bool debug,
	const std::string& debugResultsSubdir)
{
	using namespace vcl;

	if (grid.rows == 0 || grid.cols == 0 ||
		cells.size() != static_cast<size_t>(grid.rows) * grid.cols) {
		return false;
	}

	std::vector<unsigned char> exteriorWhite(cells.size(), false);
	std::vector<uint> queue;
	queue.reserve(cells.size());

	const auto addExteriorWhite = [&](uint idx) {
		if (!cells[idx].hasHit && !exteriorWhite[idx]) {
			exteriorWhite[idx] = true;
			queue.push_back(idx);
		}
	};

	//set perimeter white cells as exterior
	for (uint col = 0; col < grid.cols; ++col) {
		addExteriorWhite(col);
		addExteriorWhite((grid.rows - 1) * grid.cols + col);
	}
	for (uint row = 0; row < grid.rows; ++row) {
		addExteriorWhite(row * grid.cols);
		addExteriorWhite(row * grid.cols + grid.cols - 1);
	}

	//flood-fill from the exterior white cells to mark all connected white cells as exterior
	for (size_t readIdx = 0; readIdx < queue.size(); ++readIdx) {
		forEachCrossNeighbor(queue[readIdx], grid, [&](uint neighborIdx) {
			addExteriorWhite(neighborIdx);
			return true;
		});
	}

	//if any white cell is not marked as exterior, it is an enclosed hole
	std::vector<uint> enclosedWhiteCellIds;
	for (uint idx = 0; idx < cells.size(); ++idx) {
		if (!cells[idx].hasHit && !exteriorWhite[idx]) {
			enclosedWhiteCellIds.push_back(idx);
		}
	}

	if (enclosedWhiteCellIds.empty()) {
		return false;
	}

	if (debug) {
		PolyMesh enclosedWhiteMesh;
		enclosedWhiteMesh.enablePerVertexColor();

		for (const uint idx : enclosedWhiteCellIds) {
			const Point3d cellCenter = cells[idx].cellCenter;
			const uint vertexId = enclosedWhiteMesh.addVertex(cellCenter);
			enclosedWhiteMesh.vertex(vertexId).color() = Color::Green;
		}

		const std::filesystem::path debugOutputDir =
			std::filesystem::path(RESULTS_PATH) / debugResultsSubdir;
		std::filesystem::create_directories(debugOutputDir);
		vcl::saveMesh(
			enclosedWhiteMesh,
			(debugOutputDir / "enclosed_white_cells.ply").string());
	}

	return true;
}

#endif
