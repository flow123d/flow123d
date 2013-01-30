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
#include "mesh/region.hh"
#include "mesh/elements.h"
#include "mesh/mesh.h"

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
public:
    ElementAccessor() {}

    ElementAccessor(Mesh *mesh, unsigned int idx, bool boundary)
    : mesh_(mesh), boundary_(boundary), element_idx_(idx)
    {
       dim_=element()->dim();
    }

    inline unsigned int dim() const
        { return dim_; }

    inline const ElementIter element() const {
        if (boundary_) return &(mesh_->bc_elements[element_idx_]);
        else return &(mesh_->element[element_idx_]);
    }

    inline Region region() const
        { return element()->region(); }

    /// We need this method after replacing Region by RegionIdx, and movinf RegionDB instance into particular mesh
    inline unsigned int region_id() const {
        return element()->region_.id();
    }

    //const BoundingBox &bounding_box();

    inline bool is_boundary() const {
        return boundary_;
    }

    inline unsigned int idx() const {
        return element_idx_;
    }
private:
    /// Dimension of reference element.
    unsigned int dim_;
    // BoundingBox box_;

    /// Pointer to the mesh owning the element.
    Mesh *mesh_;
    /// True if the element is boundary, i.e. stored in Mesh::bc_elements, bulk elements are stored in Mesh::element
    bool boundary_;

    /// Index into Mesh::bc_elements or Mesh::element array.
    unsigned int element_idx_;
};




/******************************************************************* implementations
 *
 *
 */
/*
template<int spacedim>
const BoundingBox &ElementAccessor<spacedim>::bounding_box() {
    return box_;
}
*/

#endif /* ACCESSORS_HH_ */
