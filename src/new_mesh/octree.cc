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
 * $Id: octree.cc 1567 2012-02-28 13:24:58Z jan.brezina $
 * $Revision: 1567 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2012-02-28 14:24:58 +0100 (Tue, 28 Feb 2012) $
 *
 *
 */

#include "system/system.hh"
#include "new_mesh/octree.hh"
#include <typeinfo>


Octree::Octree(Mesh* _mesh)
{
	xprintf(Msg, " - Octree->Octree(Mesh* _mesh)\n");

	mesh = _mesh;
	leaf = false;
	depth = 0;

	bounding_box();
	element_boxes();
	split_area();
	distribute_elements();

	xprintf(Msg, " - Tree created\n");
}

Octree::Octree(Mesh* _mesh, arma::vec3 _minCoordinates, arma::vec3 _maxCoordinates, int _elementSize, int _splitCoor, int _depth)
{
	//xprintf(Msg, " - Octree->Octree(Mesh*, arma::vec3, arma::vec3, int, int, int)\n");

	mesh = _mesh;
	leaf = false;
	minCoordinates = _minCoordinates;
	maxCoordinates = _maxCoordinates;
	n_elements = 0;
	elements = new BoxElement * [_elementSize];
	splitCoor = _splitCoor;
	depth = _depth;
}

bool Octree::contains_element(int coor, double min, double max)
{
	return (min < maxCoordinates(coor)) & (max > minCoordinates(coor));
}

bool Octree::contains_point(arma::vec3 &_point)
{
	return IN_INTERVAL(_point(COOR_X), minCoordinates(COOR_X), maxCoordinates(COOR_X))
			& IN_INTERVAL(_point(COOR_Y), minCoordinates(COOR_Y), maxCoordinates(COOR_Y))
			& IN_INTERVAL(_point(COOR_Z), minCoordinates(COOR_Z), maxCoordinates(COOR_Z));
}

int Octree::addElement(BoxElement* _element)
{
	elements[n_elements] = _element;
	n_elements++;
	return 0;
}

int Octree::getElementCount()
{
	return n_elements;
}

void Octree::split_distribute()
{
	if (n_elements>AREA_ELEMENT_LIMIT)
	{
		split_area();
		if (!leaf) distribute_elements();
	}
	else
	{
		leaf = true;
	}
}

void Octree::bounding_box()
{
	Node* node = mesh->node_vector.begin();
	arma::vec3 point = node->point();

	minCoordinates(COOR_X) = point(COOR_X);
	minCoordinates(COOR_Y) = point(COOR_Y);
	minCoordinates(COOR_Z) = point(COOR_Z);
	maxCoordinates(COOR_X) = point(COOR_X);
	maxCoordinates(COOR_Y) = point(COOR_Y);
	maxCoordinates(COOR_Z) = point(COOR_Z);

	FOR_NODES(mesh, node )
	{
		arma::vec3 point = node->point();

		minCoordinates(COOR_X) = MIN(minCoordinates(COOR_X), point(COOR_X));
		minCoordinates(COOR_Y) = MIN(minCoordinates(COOR_Y), point(COOR_Y));
		minCoordinates(COOR_Z) = MIN(minCoordinates(COOR_Z), point(COOR_Z));
		maxCoordinates(COOR_X) = MAX(maxCoordinates(COOR_X), point(COOR_X));
		maxCoordinates(COOR_Y) = MAX(maxCoordinates(COOR_Y), point(COOR_Y));
		maxCoordinates(COOR_Z) = MAX(maxCoordinates(COOR_Z), point(COOR_Z));
	}

}

void Octree::element_boxes()
{
	n_elements = mesh->element.size();
	elements = new BoxElement * [n_elements];
	FOR_ELEMENTS(mesh, element)
	{
		arma::vec3 minCoor;
		arma::vec3 maxCoor;
		int id = element.id();
		minCoor(COOR_X) = element->node[0]->point()(COOR_X);
		minCoor(COOR_Y) = element->node[0]->point()(COOR_Y);
		minCoor(COOR_Z) = element->node[0]->point()(COOR_Z);
		maxCoor(COOR_X) = element->node[0]->point()(COOR_X);
		maxCoor(COOR_Y) = element->node[0]->point()(COOR_Y);
		maxCoor(COOR_Z) = element->node[0]->point()(COOR_Z);
		for (int i=1; i<element->n_nodes(); i++)
		{
			Node* node = element->node[i];
			minCoor(COOR_X) = MIN(minCoor(COOR_X), node->point()(COOR_X));
			minCoor(COOR_Y) = MIN(minCoor(COOR_Y), node->point()(COOR_Y));
			minCoor(COOR_Z) = MIN(minCoor(COOR_Z), node->point()(COOR_Z));
			maxCoor(COOR_X) = MAX(maxCoor(COOR_X), node->point()(COOR_X));
			maxCoor(COOR_Y) = MAX(maxCoor(COOR_Y), node->point()(COOR_Y));
			maxCoor(COOR_Z) = MAX(maxCoor(COOR_Z), node->point()(COOR_Z));
		}
		BoxElement* boxElement = new BoxElement(id, minCoor, maxCoor);
		elements[element.index()] = boxElement;
	}
}

void Octree::split_area()
{
	int medianStep = n_elements / AREA_MEDIAN_COUNT;
	int medianPosition = (int)(AREA_MEDIAN_COUNT/2);
	double median;
	double coors[AREA_MEDIAN_COUNT];
	arma::vec3 diff = maxCoordinates - minCoordinates;

	if (diff(COOR_X) >= diff(COOR_Y) && diff(COOR_X) >= diff(COOR_Z)) splitCoor = 0;
	else if (diff(COOR_Y) >= diff(COOR_Z)) splitCoor = 1;
	else splitCoor = 2;

	for (int i=0; i<AREA_MEDIAN_COUNT; i++)
	{
		coors[i] = elements[i * medianStep]->getCenterCoor(splitCoor);
	}

	for (int i=0; i<medianPosition; i++)
	{
		int minIndex=i, maxIndex=i;
		double min=coors[i], max=coors[i], change;

		for (int j=i+1; j<AREA_MEDIAN_COUNT-i; j++)
		{
			if (coors[j]<min)
			{
				min=coors[j];
				minIndex=j;
			}
			else if (coors[j]>max)
			{
				max=coors[j];
				maxIndex=j;
			}
		}
		change = coors[i];
		coors[i] = coors[minIndex];
		coors[minIndex] = change;
		if (maxIndex==i) maxIndex=minIndex;
		change = coors[AREA_MEDIAN_COUNT-i-1];
		coors[AREA_MEDIAN_COUNT-i-1] = coors[maxIndex];
		coors[maxIndex] = change;
	}

	median = coors[medianPosition];

	if (median == minCoordinates(splitCoor) || median == maxCoordinates(splitCoor))
	{
		leaf = true;
	}
	else
	{
		for (int i=0; i<CHILD_COUNT; i++)
		{
			arma::vec3 minCoor;
			arma::vec3 maxCoor;
			for (int j=0; j<3; j++)
			{
				minCoor(j) = (j==splitCoor && i==1) ? median : minCoordinates(j);
				maxCoor(j) = (j==splitCoor && i==0) ? median : maxCoordinates(j);
			}

			child[i] = new Octree(mesh, minCoor, maxCoor, n_elements, splitCoor, depth+1);
		}
	}
}

void Octree::distribute_elements()
{
	for (int i=0; i<n_elements; i++)
	{
		for (int j=0; j<CHILD_COUNT; j++)
		{
			if (child[j]->contains_element(splitCoor, elements[i]->getMinCoor(splitCoor), elements[i]->getMaxCoor(splitCoor)))
			{
				child[j]->addElement(elements[i]);
			}
		}
	}
	child[0]->split_distribute();
	child[1]->split_distribute();
}
