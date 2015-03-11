/*
 * bc_field.hh
 *
 *  Created on: Feb 13, 2014
 *      Author: jb
 */

#ifndef BC_FIELD_HH_
#define BC_FIELD_HH_


#include "field.hh"


/**
 * Same as Field<...> but for boundary regions.
 *
 * Definition of BCField must be in separate file.
 * In other case source file field.cc is too big and compiler can throw compile error.
 */
template<int spacedim, class Value>
class BCField : public Field<spacedim, Value> {
public:
    BCField() : Field<spacedim,Value>("anonymous_bc", true) {}
};

#endif /* BC_FIELD_HH_ */
