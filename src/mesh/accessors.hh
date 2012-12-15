/*
 * Accessors.hh
 *
 *  Created on: Dec 4, 2012
 *      Author: jb
 */

#ifndef ACCESSORS_HH_
#define ACCESSORS_HH_

#include "new_mesh/bounding_box.hh"
#include "mesh/mesh_types.hh"

/**
 * Element accessor templated just by dimension of the embedding space, used by Fields.
 * This should allow algorithms over elements where dimension of particular element is runtime parameter.
 *
 * TODO: (add various things needed by particular Field implementations)
 *
 * need function to calculate intersection (object) of two ElementAccessors, but this definitely should be templated by
 * dimension of the ref. element (or rather shape of ref. element), here we can have case dispatch
 *
 */
template <int spacedim>
class ElementAccessor {
    inline unsigned int dim() const
        { return dim_; }
    inline const ElementIter element() const
        { return element_;}

    const BoundingBox &bounding_box();

private:
    /// Dimension of reference element.
    unsigned int dim_;
    BoundingBox box_;

    // TODO: remove pointer reference as soon as possible (provide other necessary access functions for FieldInterpoleted P0)
    ElementIter element_;
};




/******************************************************************* implementations
 *
 *
 */

template<int spacedim>
const BoundingBox &ElementAccessor<spacedim>::bounding_box() {
    return box_;
}


#endif /* ACCESSORS_HH_ */
