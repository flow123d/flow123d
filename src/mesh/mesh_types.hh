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

class Node;
class Element;
class Boundary;
class Edge;

// Preparation for next development
typedef flow::VectorId<Node> NodeVector;
typedef NodeVector::Iter NodeIter;
typedef NodeVector::FullIter NodeFullIter;

// iterator over elements
// should be mesh member, but then we have problem how to have ElementIter as memeber of
// Node or other classes without cyclic inclusion
typedef flow::VectorId<Element> ElementVector;
typedef ElementVector::Iter ElementIter;
typedef ElementVector::FullIter ElementFullIter;

typedef flow::Vector<Boundary> BoundaryVector;
typedef BoundaryVector::Iter BoundaryIter;
typedef BoundaryVector::FullIter BoundaryFullIter;

typedef flow::Vector<Edge> EdgeVector;
typedef EdgeVector::Iter EdgeIter;
typedef EdgeVector::FullIter EdgeFullIter;

#endif /* MESH_TYPES_HH_ */
