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
 * $Id: octree.hh 1567 2012-02-28 13:24:58Z jan.brezina $
 * $Revision: 1567 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2012-02-28 14:24:58 +0100 (Tue, 28 Feb 2012) $
 *
 *
 */

#ifndef OCTREE_HH_
#define OCTREE_HH_

#include "system/system.hh"
#include "mesh/mesh.h"
#include "new_mesh/box_element.hh"
#include <armadillo>

#define MIN(x, y) (x < y ? x : y)
#define MAX(x, y) (x > y ? x : y)
#define IN_INTERVAL(x, down, up) ((x >= down) & (x <= up) ? true : false)
#define AREA_ELEMENT_LIMIT 20
#define CHILD_COUNT 2 //don't change
#define AREA_MEDIAN_COUNT 9 //must be even
#define COOR_X 0
#define COOR_Y 1
#define COOR_Z 2

class Octree {
public:
    Octree(Mesh* _mesh);
    Octree(Mesh* _mesh, arma::vec3 _minCoordinates, arma::vec3 _maxCoordinates, int _elementSize, int _splitCoor, int _depth);
    bool contains_point(arma::vec3 &_point);
    bool contains_element(int coor, double min, double max);
    int addElement(BoxElement* _element);
    int getElementCount();
    void split_distribute();

private:
    void bounding_box();
    void element_boxes();
    void split_area();
    void distribute_elements();
    Mesh* mesh;
    Octree* child[CHILD_COUNT];
    BoxElement** elements;
    int n_elements;
    int splitCoor;
    int depth;
    bool leaf;
    arma::vec3 minCoordinates;
    arma::vec3 maxCoordinates;
};


#endif /* OCTREE_HH_ */
