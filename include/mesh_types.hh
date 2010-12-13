/*
 * mesh_types.hh
 *
 *  Created on: Aug 6, 2010
 *      Author: jb
 */

#ifndef MESH_TYPES_HH_
#define MESH_TYPES_HH_

#include <sys_vector.hh>

class Node;
class Element;
class Boundary;

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

typedef flow::VectorId<Boundary> BoundaryVector;
typedef BoundaryVector::Iter BoundaryIter;
typedef BoundaryVector::FullIter BoundaryFullIter;

#endif /* MESH_TYPES_HH_ */
