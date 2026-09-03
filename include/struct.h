#ifndef VCL_TEST_EXTERNAL_888_MOLD_CHECK_STRUCT_H
#define VCL_TEST_EXTERNAL_888_MOLD_CHECK_STRUCT_H

#include <array>
#include <vector>

#include <vclib/space/core/point.h>

struct GridChoice
{
	vcl::uint rows = 1;
	vcl::uint cols = 1;

	double sideU = 0.0;
	double sideV = 0.0;

	double minU = 0.0;
	double minV = 0.0;
	double maxU = 0.0;
	double maxV = 0.0;
};

struct CellData
{
	std::array<vcl::Point3d, 4> cellCorners;
	vcl::Point3d cellCenter;
	double distance = 0.0; //distance from plane
	std::array<double, 2> boundaries = {0.0, 0.0}; // min and max distance from mold plane
	std::vector<vcl::Point3d> hitPoints; //all mesh hitpoints on the cell, if any
	bool hasHit = false; // true if the cell has at least one hit point
	bool hasClampedHit = false; // true if cell does not respect draft angle with hit point
	bool isReduced = false; // true if is removed from ReducePoints
	bool isBiharmonicFilledHit = false; // true if the cell is filled by biharmonic interpolation
	bool isMovedForward = false; // cyan points
	bool isDiscarded = false; // points outside box constraints
	bool isInside = false; // inside the mesh, used for biharmonic interpolation
	// double clampedDistance = 0.0;
};

struct HitCellShapeData
{
	double area = 0.0;
	double perimeter = 0.0;
	double compactness = 0.0;
};

#endif
