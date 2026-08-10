#ifndef VCL_TEST_EXTERNAL_888_MOLD_CHECK_RAYS_H
#define VCL_TEST_EXTERNAL_888_MOLD_CHECK_RAYS_H

#include "helper.h"
#include "struct.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <unordered_set>
#include <vector>

#include <vclib/embree/scene.h>
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

#endif