#ifndef VCL_TEST_EXTERNAL_888_MOLD_CHECK_DEBUG_OUTPUT_H
#define VCL_TEST_EXTERNAL_888_MOLD_CHECK_DEBUG_OUTPUT_H

#include "struct.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <filesystem>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

#include <vclib/io.h>
#include <vclib/meshes.h>

static void validateClampedCells(
	const std::vector<CellData>& clampedCells,
	const std::vector<vcl::uint>& allCells,
	const vcl::Point3d& direction,
	double coneCosThreshold,
	float eps)
{
	using namespace vcl;

	std::atomic<uint> violatingPoints{0};

	vcl::parallelFor(allCells, [&](uint i) {
		if (!clampedCells[i].hasHit) {
			return;
		}

		const Point3d& point = clampedCells[i].hitPoints[0];

		for (uint j = 0; j < clampedCells.size(); ++j) {
			if (i == j || !clampedCells[j].hasHit) {
				continue;
			}

			const Point3d dirToOther =
				clampedCells[j].hitPoints[0] - point;

			const double norm = dirToOther.norm();

			if (norm < eps) {
				continue;
			}

			const double cosVal =
				(dirToOther / norm).dot(-direction);

			if (cosVal > coneCosThreshold + eps) {
				violatingPoints.fetch_add(1);
				return;
			}
		}
	});

	const uint totalViolations = violatingPoints.load();

	std::cout << "Clamped validation: "
			  << (totalViolations == 0 ? "OK" : "VIOLATIONS")
			  << " (violating points: " << totalViolations << ")\n";

	std::cout.flush();
}

static vcl::TriMesh createMoldSurface(
	const std::vector<CellData>& clampedCells,
	const GridChoice& grid,
	const vcl::Point3d& direction)
{
	using namespace vcl;

	TriMesh tm;
	tm.enablePerFaceColor();

	for (uint row = 0; row + 1 < grid.rows; row += 1) {
		for (uint col = 0; col + 1 < grid.cols; col += 1) {
			const uint c00 = row * grid.cols + col;
			const uint c10 = c00 + 1;
			const uint c01 = c00 + grid.cols;
			const uint c11 = c01 + 1;

			const std::array<const CellData*, 4> cells = {
				&clampedCells[c00],
				&clampedCells[c10],
				&clampedCells[c01],
				&clampedCells[c11]};

			const double averageDistance =
				(cells[0]->distance +
				 cells[1]->distance +
				 cells[2]->distance +
				 cells[3]->distance) *
				0.25;

			const Point3d foot =
				(cells[0]->cellCenter +
				 cells[1]->cellCenter +
				 cells[2]->cellCenter +
				 cells[3]->cellCenter) *
				0.25;

			const Point3d medianPoint = foot + direction * averageDistance;

			const Point3d p0 = cells[0]->hitPoints[0];
			const Point3d p1 = cells[1]->hitPoints[0];
			const Point3d p2 = cells[2]->hitPoints[0];
			const Point3d p3 = cells[3]->hitPoints[0];

			const uint v0 = tm.addVertex(p0);
			const uint v1 = tm.addVertex(p1);
			const uint v2 = tm.addVertex(p2);
			const uint v3 = tm.addVertex(p3);
			const uint vc = tm.addVertex(medianPoint);

			tm.addFace(v0, vc, v1);
			tm.addFace(v1, vc, v3);
			tm.addFace(v3, vc, v2);
			tm.addFace(v2, vc, v0);
		}
	}

	return tm;
}

static vcl::uint addFaceWithColor(
	vcl::TriMesh& tm,
	vcl::uint v0,
	vcl::uint v1,
	vcl::uint v2,
	const vcl::Color& faceColor)
{
	tm.enablePerFaceColor();
	const vcl::uint fid = tm.addFace(v0, v1, v2);
	tm.face(fid).color() = faceColor;
	return fid;
}

static void addQuadPrism(
	vcl::TriMesh& tm,
	const std::array<vcl::Point3d, 4>& baseCorners,
	double startOffset,
	double endOffset,
	const vcl::Point3d& dir,
	const vcl::Color& faceColor)
{
	using namespace vcl;

	tm.enablePerFaceColor();

	std::array<Point3d, 4> b;
	std::array<Point3d, 4> t;

	for (uint k = 0; k < 4; ++k) {
		b[k] = baseCorners[k] + dir * startOffset;
		t[k] = baseCorners[k] + dir * endOffset;
	}

	std::array<uint, 8> ids;

	for (uint k = 0; k < 4; ++k) {
		ids[k + 0] = tm.addVertex(b[k]);
		ids[k + 4] = tm.addVertex(t[k]);
	}

	addFaceWithColor(tm, ids[0], ids[2], ids[1], faceColor);
	addFaceWithColor(tm, ids[0], ids[3], ids[2], faceColor);

	addFaceWithColor(tm, ids[4], ids[5], ids[6], faceColor);
	addFaceWithColor(tm, ids[4], ids[6], ids[7], faceColor);

	addFaceWithColor(tm, ids[0], ids[1], ids[5], faceColor);
	addFaceWithColor(tm, ids[0], ids[5], ids[4], faceColor);

	addFaceWithColor(tm, ids[1], ids[2], ids[6], faceColor);
	addFaceWithColor(tm, ids[1], ids[6], ids[5], faceColor);

	addFaceWithColor(tm, ids[2], ids[3], ids[7], faceColor);
	addFaceWithColor(tm, ids[2], ids[7], ids[6], faceColor);

	addFaceWithColor(tm, ids[3], ids[0], ids[4], faceColor);
	addFaceWithColor(tm, ids[3], ids[4], ids[7], faceColor);
}


static double pointOffsetFromCellPlane(
	const CellData& cell,
	const vcl::Point3d& point,
	const vcl::Point3d& direction)
{
	return (point - cell.cellCenter).dot(direction);
}

static vcl::TriMesh createRemainingMold(
	const std::vector<CellData>& cells,
	const std::vector<CellData>& clampedCells,
	const std::vector<CellData>& moldClampedCells,
	const vcl::Point3d& direction)
{
	using namespace vcl;

	TriMesh tm;

	if (moldClampedCells.size() != clampedCells.size()) {
		return tm;
	}

	for (uint i = 0; i < clampedCells.size(); ++i) {
		if (!clampedCells[i].hasHit || !moldClampedCells[i].hasHit) {
			continue;
		}

		const double moldDistance =
			pointOffsetFromCellPlane(
				clampedCells[i],
				moldClampedCells[i].hitPoints[0],
				direction);

		if (cells[i].hasHit) {
			for (std::size_t hitIndex = 2;
				 hitIndex + 1 < cells[i].hitPoints.size() - 1;
				 hitIndex += 2) {
				const double startHitDistance =
					pointOffsetFromCellPlane(
						clampedCells[i],
						cells[i].hitPoints[hitIndex],
						direction);
				const double endHitDistance =
					pointOffsetFromCellPlane(
						clampedCells[i],
						cells[i].hitPoints[hitIndex + 1],
						direction);
				addQuadPrism(
					tm,
					clampedCells[i].cellCorners,
					startHitDistance,
					endHitDistance,
					direction,
					Color::Red);
			}

			const double lastHitDistance =
				pointOffsetFromCellPlane(
					clampedCells[i],
					clampedCells[i].hitPoints.back(),
					direction);
			addQuadPrism(
				tm,
				clampedCells[i].cellCorners,
				lastHitDistance,
				moldDistance,
				direction,
				Color::Red);

			if (cells[i].distance != clampedCells[i].distance) {
				addQuadPrism(
					tm,
					clampedCells[i].cellCorners,
					clampedCells[i].distance,
					cells[i].distance,
					direction,
					Color::Red);
			}
		}
		else {
			addQuadPrism(
				tm,
				clampedCells[i].cellCorners,
				clampedCells[i].distance,
				moldDistance,
				direction,
				Color::Red);
		}
	}

	return tm;
}

static void addSegment(
	vcl::EdgeMesh& em,
	const vcl::Point3d& a,
	const vcl::Point3d& b)
{
	const vcl::uint va = em.addVertex(a);
	const vcl::uint vb = em.addVertex(b);

	em.addEdge(va, vb);
}

static void addColoredPoint(
	vcl::PolyMesh& mesh,
	const vcl::Point3d& point,
	const vcl::Color& color)
{
	const vcl::uint vId = mesh.addVertex(point);
	mesh.vertex(vId).color() = color;
}

static vcl::TriMesh makeDebugPlaneMesh(
	const GridChoice& grid,
	const vcl::Point3d& planePoint,
	const vcl::Point3d& u,
	const vcl::Point3d& v)
{
	using namespace vcl;

	TriMesh planeMesh;

	const Point3d p0 =
		planePoint + u * grid.minU + v * grid.minV;

	const Point3d p1 =
		planePoint + u * grid.maxU + v * grid.minV;

	const Point3d p2 =
		planePoint + u * grid.maxU + v * grid.maxV;

	const Point3d p3 =
		planePoint + u * grid.minU + v * grid.maxV;

	const uint v0 = planeMesh.addVertex(p0);
	const uint v1 = planeMesh.addVertex(p1);
	const uint v2 = planeMesh.addVertex(p2);
	const uint v3 = planeMesh.addVertex(p3);

	planeMesh.addFace(v0, v1, v2);
	planeMesh.addFace(v0, v2, v3);

	return planeMesh;
}

static vcl::EdgeMesh createPerimeterSegments(
	const std::vector<vcl::uint>& componentIndices,
	const std::vector<CellData>& cells,
	const GridChoice& grid)
{
	using namespace vcl;

	EdgeMesh em;

	std::unordered_set<uint> componentSet(
		componentIndices.begin(), componentIndices.end());

	for (uint idx : componentIndices) {
		const uint row = idx / grid.cols;
		const uint col = idx % grid.cols;
		const CellData& cell = cells[idx];

		if (col == 0 || componentSet.count(idx - 1) == 0) {
			addSegment(em, cell.cellCorners[0], cell.cellCorners[3]);
		}
		if (col + 1 == grid.cols || componentSet.count(idx + 1) == 0) {
			addSegment(em, cell.cellCorners[1], cell.cellCorners[2]);
		}
		if (row == 0 || componentSet.count(idx - grid.cols) == 0) {
			addSegment(em, cell.cellCorners[0], cell.cellCorners[1]);
		}
		if (row + 1 == grid.rows || componentSet.count(idx + grid.cols) == 0) {
			addSegment(em, cell.cellCorners[3], cell.cellCorners[2]);
		}
	}

	return em;
}

static void debugOutput(
	const std::vector<CellData>& cells,
	const std::vector<CellData>& clampedCells,
	const std::vector<CellData>& moldClampedCells,
	const ConnectedComponentData& largestComponent,
	const GridChoice& grid,
	const vcl::Point3d& planePoint,
	const vcl::Point3d& u,
	const vcl::Point3d& v,
	const vcl::Point3d& direction,
	float eps,
	const std::string& debugResultsSubdir,
	double totalAreaHit,
	double clampedAreaHit,
	double percentClamped,
	double hiddenAreaHit,
	double percentHidden,
	double componentRatio,
	double qualityScore)
{
	using namespace vcl;

	PolyMesh hitPointsMesh;
	hitPointsMesh.enablePerVertexColor();
	for (uint i = 0; i < cells.size(); ++i) {
		if (!cells[i].hasHit) continue;
		addColoredPoint(hitPointsMesh, cells[i].hitPoints[0], Color::Yellow);
	}
	
	PolyMesh clampedonlyPointsMesh;
	clampedonlyPointsMesh.enablePerVertexColor();
	for (uint i = 0; i < clampedCells.size(); ++i) {
		if (!cells[i].hasHit) continue;
		if (cells[i].distance == clampedCells[i].distance) continue;
		addColoredPoint(clampedonlyPointsMesh, clampedCells[i].hitPoints[0], Color::Red);
	}
	
	PolyMesh clampednohitPointsMesh;
	clampednohitPointsMesh.enablePerVertexColor();
	for (uint i = 0; i < clampedCells.size(); ++i) {
		if (cells[i].hasHit) continue;
		if (cells[i].distance == clampedCells[i].distance) continue;
		addColoredPoint(clampednohitPointsMesh, clampedCells[i].hitPoints[0], Color::White);
	}
	
	PolyMesh clampedPointsMesh;
	clampedPointsMesh.enablePerVertexColor();
	for (uint i = 0; i < clampedCells.size(); ++i) {
		addColoredPoint(clampedPointsMesh, clampedCells[i].hitPoints[0], Color::Blue);
	}

	PolyMesh moldClampedPointsMesh;
	moldClampedPointsMesh.enablePerVertexColor();
	for (uint i = 0; i < moldClampedCells.size(); ++i) {
		if (!moldClampedCells[i].hasHit) continue;
		addColoredPoint(moldClampedPointsMesh, moldClampedCells[i].hitPoints[0], Color::Cyan);
	}

	PolyMesh thirdHitPointsMesh;
	thirdHitPointsMesh.enablePerVertexColor();
	for (uint i = 0; i < cells.size(); ++i) {
		if (cells[i].hitPoints.size() <= 2) continue;
		addColoredPoint(thirdHitPointsMesh, cells[i].hitPoints[2], Color::Magenta);
	}

	PolyMesh lastHitPointsMesh;
	lastHitPointsMesh.enablePerVertexColor();
	for (uint i = 0; i < cells.size(); ++i) {
		if (!cells[i].hasHit) continue;
		addColoredPoint(lastHitPointsMesh, cells[i].hitPoints.back(), Color::Green);
	}

	PolyMesh oddHitPointsMesh;
	oddHitPointsMesh.enablePerVertexColor();
	uint oddHitCells = 0;
	for (uint i = 0; i < cells.size(); ++i) {
		if (!cells[i].hasHit) continue;
		if (cells[i].hitPoints.size() % 2 == 0) continue;
		++oddHitCells;
		for (const Point3d& hitPoint : cells[i].hitPoints) {
			addColoredPoint(oddHitPointsMesh, hitPoint, Color::LightYellow);
		}
	}
	
	PolyMesh missedPointsMesh;
	missedPointsMesh.enablePerVertexColor();
	for (uint i = 0; i < clampedCells.size(); ++i) {
		if (clampedCells[i].hasHit) continue;
		addColoredPoint(missedPointsMesh, clampedCells[i].hitPoints[0], Color::Green);
	}

	PolyMesh largestComponentMesh;
	largestComponentMesh.enablePerVertexColor();
	for (uint i : largestComponent.indices) {
		addColoredPoint(largestComponentMesh, clampedCells[i].hitPoints[0], Color::Cyan);
	}

	const TriMesh planeMesh =
		makeDebugPlaneMesh(grid, planePoint, u, v);

	TriMesh ClampedPrismMesh;
	for (uint i = 0; i < clampedCells.size(); ++i) {
		if (!clampedCells[i].hasHit) continue;
		addQuadPrism(ClampedPrismMesh, clampedCells[i].cellCorners, -eps, clampedCells[i].distance, direction, vcl::Color::White);
	}
	
	TriMesh remainingMoldMesh;
	EdgeMesh segmentsRemainingMold;
	for (uint i = 0; i < clampedCells.size(); ++i) {
		if (!cells[i].hasHit) continue;
		if (cells[i].distance == clampedCells[i].distance) continue;
		addQuadPrism(remainingMoldMesh, clampedCells[i].cellCorners, clampedCells[i].distance, cells[i].distance, direction, vcl::Color::Red);
		addSegment(
			segmentsRemainingMold,
			clampedCells[i].hitPoints[0],
			cells[i].hitPoints[0]);
	}

	const TriMesh moldSurfaceMesh = createMoldSurface(clampedCells, grid, direction);

	const EdgeMesh perimeterSegments =
		createPerimeterSegments(
			largestComponent.indices, clampedCells, grid);

	
	const std::filesystem::path debugOutputDir =
		std::filesystem::path(RESULTS_PATH) /
		debugResultsSubdir;

	std::filesystem::create_directories(debugOutputDir);

	for (const auto& entry :
		 std::filesystem::directory_iterator(debugOutputDir)) {
		if (entry.is_regular_file() && entry.path().extension() == ".ply") {
			std::filesystem::remove(entry.path());
		}
	}

	const std::string base =
		(debugOutputDir / "888_mold_check").string();
	saveMesh(hitPointsMesh, base + "_hit_points.ply");
	saveMesh(clampedonlyPointsMesh, base + "_clamped_only_points.ply");
	saveMesh(clampednohitPointsMesh, base + "_clamped_nohit_points.ply");
	saveMesh(clampedPointsMesh, base + "_all_clamped_points.ply");
	saveMesh(moldClampedPointsMesh, base + "_mold_clamped_points.ply");
	saveMesh(thirdHitPointsMesh, base + "_third_hit_points.ply");
	saveMesh(lastHitPointsMesh, base + "_last_hit_points.ply");
	saveMesh(oddHitPointsMesh, base + "_odd_hit_points.ply");
	saveMesh(planeMesh, base + "_plane.ply");
	saveMesh(missedPointsMesh, base + "_missed_points.ply");
	saveMesh(ClampedPrismMesh, base + "_clamped_prisms.ply");
	saveMesh(remainingMoldMesh, base + "_remaining_mold.ply");
	saveMesh(segmentsRemainingMold, base + "_remaining_mold_segments.ply");
	saveMesh(moldSurfaceMesh, base + "_mold_surface.ply");
	saveMesh(largestComponentMesh, base + "_largest_component_points.ply");
	saveMesh(perimeterSegments, base + "_largest_component_perimeter.ply");
	
	std::cout << "Clamped points: " << clampedPointsMesh.vertexCount() << "\n";
	std::cout << "Mold clamped points: " << moldClampedPointsMesh.vertexCount() << "\n";
	std::cout << "Third hit points: " << thirdHitPointsMesh.vertexCount() << "\n";
	std::cout << "Odd hit cells: " << oddHitCells << "\n";
	std::cout << "Odd hit points: " << oddHitPointsMesh.vertexCount() << "\n";
	std::cout << "Mold surface median points: " << moldSurfaceMesh.vertexCount() << "\n";
	std::cout << "Largest component cells: "
			  << largestComponent.indices.size() << "\n";
	std::cout << "Largest component area: "
			  << largestComponent.area << "\n";
	std::cout << "Largest component perimeter: "
			  << largestComponent.perimeter << "\n";
	std::cout << "Largest component compactness: "
			  << largestComponent.compactness << "\n";
	std::cout << "TotalAreaHit: "
			  << totalAreaHit << "\n";
	std::cout << "ClampedAreaHit: "
			  << clampedAreaHit << "\n";
	std::cout << "percentClamped: "
			  << percentClamped << "\n";
	std::cout << "hiddenAreaHit: "
			  << hiddenAreaHit << "\n";
	std::cout << "percentHidden: "
			  << percentHidden << "\n";
	std::cout << "componentRatio: "
			  << componentRatio << "\n";
	std::cout << "qualityScore: "
			  << qualityScore << "\n";
	std::cout << "Saved debug meshes:\n"
			<< " - " << base << "_hit_points.ply\n"
			<< " - " << base << "_clamped_only_points.ply\n"
			<< " - " << base << "_clamped_nohit_points.ply\n"
			<< " - " << base << "_all_clamped_points.ply\n"
			<< " - " << base << "_mold_clamped_points.ply\n"
			<< " - " << base << "_third_hit_points.ply\n"
			<< " - " << base << "_last_hit_points.ply\n"
			<< " - " << base << "_odd_hit_points.ply\n"
			<< " - " << base << "_plane.ply\n"
			<< " - " << base << "_missed_points.ply\n"
			<< " - " << base << "_clamped_prisms.ply\n"
			<< " - " << base << "_remaining_mold.ply\n"
			<< " - " << base << "_remaining_mold_segments.ply\n"
			<< " - " << base << "_mold_surface.ply\n"
			<< " - " << base << "_largest_component_points.ply\n"
			<< " - " << base << "_largest_component_perimeter.ply\n";

	std::cout << "=== moldCheck completed successfully ===\n";
	std::cout.flush();
}

#endif
