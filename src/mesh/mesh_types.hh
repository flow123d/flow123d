/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    mesh_types.hh
 * @brief   
 * @author  Jan Brezina
 */

#ifndef MESH_TYPES_HH_
#define MESH_TYPES_HH_

#include "system/sys_vector.hh"
#include <vector>

// Forward declarations
class Node;
class Element;
class Edge;

// Preparation for next development
typedef flow::VectorId<Node> NodeVector;
typedef NodeVector::Iter NodeIter;
typedef NodeVector::FullIter NodeFullIter;

// iterator over elements
// should be mesh member, but then we have problem how to have ElementIter as memeber of
// Node or other classes without cyclic inclusion
// temporary iterator for access to vector of all (bulk + boundary) elements
typedef vector<Element>::iterator ElementIterator;

typedef flow::Vector<Edge> EdgeVector;

#endif /* MESH_TYPES_HH_ */
