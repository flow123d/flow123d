/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id: interpolation_main.cc 1567 2012-02-28 13:24:58Z jan.brezina $
 * $Revision: 1567 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2012-02-28 14:24:58 +0100 (Tue, 28 Feb 2012) $
 *
 *
 */

//#include "system/sys_profiler.hh"
#include "system/system.hh"
#include "mesh/msh_gmshreader.h"
//#include "new_mesh/mesh.hh"
#include "mesh/mesh.h"
#include "new_mesh/bounding_interval_hierarchy.hh"
#include "new_mesh/ngh/include/triangle.h"
#include "new_mesh/ngh/include/tetrahedron.h"
#include "new_mesh/ngh/include/polygon.h"
#include "new_mesh/ngh/include/intersection.h"
#include <armadillo>


int main(int argc, char **argv) {

	const std::string& file_name = "../tests/11_reactional_transport_semchem/input/sit_trans.msh";

	MeshReader* meshReader = new GmshMeshReader();

	const std::string& ngh_fname = "../tests/11_reactional_transport_semchem/input/sit_trans.ngh";
	const std::string& bcd_fname = "../tests/11_reactional_transport_semchem/input/sit_trans.fbc";
	Mesh* mesh = new Mesh(ngh_fname, bcd_fname);
	meshReader->read(file_name, mesh);
	BoundingIntevalHierachy* bihTree = new BoundingIntevalHierachy(mesh);

	TPoint* pointA = new TPoint(-0.10, -0.10, 0.00);
	TPoint* pointB = new TPoint(1.60, 1.60, 0.00);
	TPoint* pointC = new TPoint(0.10, 0.10, 0.50);
	TTriangle triangle(pointA, pointB, pointC);

	TPoint* point0 = new TPoint(0.00, 0.00, 0.00);
	TPoint* point1 = new TPoint(3.00, 0.00, 0.00);
	TPoint* point2 = new TPoint(0.00, 3.00, 0.00);
	TPoint* point3 = new TPoint(0.00, 0.00, 3.00);
	TTetrahedron tetrahedron(point0, point1, point2, point3);

	TIntersectionType it;
	double area;
	GetIntersection(triangle, tetrahedron, it, area);
	xprintf(Msg, "Area: %f, Triangle: %f\n", area, triangle.GetArea());

	//std::vector<BoundingBox *> searchedElements;
	//bihTree->find_elements(*triangle, searchedElements);

	/*xprintf(Msg, " - searched elements ids: ");

	for (std::vector<BoundingBox *>::iterator tmp = searchedElements.begin(); tmp!=searchedElements.end(); tmp++)
	{
		BoundingBox* b = *tmp;

		xprintf(Msg, "%d ", b->getId());
	}
	xprintf(Msg, "\n");*/

	/*TPoint* point1 = new TPoint( 1.0, 1.0, 1.0);
	TPoint* point2 = new TPoint( 2.0, 6.0, 3.0);
	TPoint* point3 = new TPoint( 2.0, 2.0, 1.0);
	TPoint* point4 = new TPoint(-1.0, 3.0, 3.0);
	TPoint* point5 = new TPoint(-1.0, 5.0, 4.0);
	TPoint* point6 = new TPoint( 3.0, 4.0, 1.5);
	TPoint* point7 = new TPoint( 0.0, 7.0, 4.5);
	TPolygon* pol = new TPolygon();
	pol->Add(point1);
	pol->Add(point2);
	pol->Add(point3);
	pol->Add(point4);
	pol->Add(point5);
	pol->Add(point6);
	pol->Add(point7);
	pol->Write();
	xprintf(Msg, "Polygon: %f\n", pol->GetArea());*/

	xprintf(Msg, " - interpolation_main executed\n");

	getchar();
}
