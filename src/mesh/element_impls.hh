/*
 * element_impls.hh
 *
 *  Created on: Jun 27, 2012
 *      Author: jb
 */

#ifndef ELEMENT_IMPLS_HH_
#define ELEMENT_IMPLS_HH_

#include "sides.h"
#include "side_impl.hh"
//#include "mesh/accessors.hh"

//template <int spacedim>
//class ElementAccessor;

inline unsigned int Element::dim() const {
    return dim_;
}


inline unsigned int Element::index() const {
    return mesh_->element.index( this );
}


inline unsigned int Element::n_nodes() const {
    return dim()+1;
}



inline unsigned int Element::n_sides() const {
    return dim()+1;
}



inline SideIter Element::side(const unsigned int loc_index) {
    return SideIter( Side(this, loc_index) );
}

/*
inline ElementAccessor<3> Element::element_accessor(){
  return mesh_->element_accessor( mesh_->element.index(this) );
}
//*/
#endif /* ELEMENT_IMPLS_HH_ */
