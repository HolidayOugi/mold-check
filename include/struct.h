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
	double distance;
	std::vector<vcl::Point3d> hitPoints;
	bool hasHit = false;
};

struct ConnectedComponentData
{
	std::vector<vcl::uint> indices;
	double area = 0.0;
	double perimeter = 0.0;
	double compactness = 0.0;
};

#endif
