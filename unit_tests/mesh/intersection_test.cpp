#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include "system/system.hh"
#include "mesh/ngh/include/point.h"
#include "mesh/ngh/include/intersection.h"


TEST(intersection, triangle_tetrahedron) {
	TIntersectionType it;
	double area;

	// create tetrahedron
	TPoint point0(0.00, 0.00, 0.00);
	TPoint point1(3.00, 0.00, 0.00);
	TPoint point2(0.00, 3.00, 0.00);
	TPoint point3(0.00, 0.00, 3.00);
	TTetrahedron tetrahedron(point0, point1, point2, point3);

	// triangle is in tetrahedron
	TPoint pointA(0.50, 0.50, 0.50);
	TPoint pointB(0.50, 1.50, 0.50);
	TPoint pointC(1.50, 0.50, 0.50);
	TTriangle triangle(pointA, pointB, pointC);

	MessageOut() << "Test - triangle is in tetrahedron\n";
	GetIntersection(triangle, tetrahedron, it, area);
	EXPECT_FLOAT_EQ(area, 0.5);

	// triangle is greater than tetrahedron, intersection is triangle
	pointA.SetCoord(-3.0, 2.0, 2.0);
	pointB.SetCoord(2.0, -3.0, 2.0);
	pointC.SetCoord(2.0, 2.0, 2.0);
	triangle.SetPoints(pointA, pointB, pointC);

	MessageOut() << "Test - triangle is greater than tetrahedron, intersection is triangle\n";
	GetIntersection(triangle, tetrahedron, it, area);
	EXPECT_FLOAT_EQ(area, 0.5);

	// intersection is tetragon
	pointA.SetCoord(-0.50, 0.50, 1.00);
	pointB.SetCoord(2.00, 0.50, 1.00);
	pointC.SetCoord(2.00, 3.00, 1.00);
	triangle.SetPoints(pointA, pointB, pointC);

	MessageOut() << "Test - intersection is tetragon\n";
	GetIntersection(triangle, tetrahedron, it, area);
	EXPECT_FLOAT_EQ(area, 0.875);

	// intersection is pentagon
	pointA.SetCoord(-1.0, 2.00, 1.00);
	pointB.SetCoord(1.50, 2.00, 1.00);
	pointC.SetCoord(1.50,-0.5 , 1.00);
	triangle.SetPoints(pointA, pointB, pointC);

	MessageOut() << "Test - intersection is pentagon\n";
	GetIntersection(triangle, tetrahedron, it, area);
	EXPECT_FLOAT_EQ(area, 0.5+0.5+0.5*3.0/4.0);

	// intersection is hexagon (plane parallel to x-y)
	pointA.SetCoord(2.00, 2.00, 0.50);
	pointB.SetCoord(2.00, -1.00, 0.50);
	pointC.SetCoord(-1.00, 2.00, 0.50);
	triangle.SetPoints(pointA, pointB, pointC);

	MessageOut() << "Test - intersection is hexagon (plane parallel to x-y)\n";
	GetIntersection(triangle, tetrahedron, it, area);
	EXPECT_FLOAT_EQ(area, 2.375);

	// intersection is hexagon
	pointA.SetCoord(0.25, 2.00, 1.00);
	pointB.SetCoord(2.25, 1.00, -1.00);
	pointC.SetCoord(0.25, -1.00, 3.00);
	triangle.SetPoints(pointA, pointB, pointC);

	MessageOut() << "Test - intersection is hexagon\n";
	GetIntersection(triangle, tetrahedron, it, area);
	EXPECT_FLOAT_EQ(area, 3.477919);
}
