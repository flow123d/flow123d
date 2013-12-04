/*
 * field_base.cc
 *
 *  Created on: Dec 4, 2012
 *      Author: jb
 */



#include "fields/field_base_impl.hh"
#include <limits>


/****************************************************************************
 *  Implementation of FieldCommon
 */

FieldCommonBase::FieldCommonBase(bool bc)
: changed_during_set_time(false),       // reading this variable is in fact invalid up to the first call of the set_time
  bc_(bc),
  n_comp_(0),
  element_selection_(NULL),
  default_( IT::Default::obligatory()),
  mesh_(NULL),
  changed_from_last_set_time_(false),
  last_set_time_( -numeric_limits<double>::infinity() )
{}




// getters
const std::string & FieldCommonBase::name() const
{ return name_; }


const std::string & FieldCommonBase::units() const
{ return units_; }


double FieldCommonBase::time() const
{ return last_set_time_; }


const std::string  FieldCommonBase::desc() const
{ if(default_.has_value_at_declaration())
    return "Default Field value: " + default_.value() + " \n " + desc_; 
  else
    return desc_;
}


const IT::Default & FieldCommonBase::get_default() const
{ return default_; }


bool FieldCommonBase::is_bc() const
{ return bc_; }


bool FieldCommonBase::is_enum_valued() const
{ return enum_valued_; }


unsigned int FieldCommonBase::n_comp() const
{ return n_comp_; }


Mesh * FieldCommonBase::mesh() const
{ return mesh_; }


bool FieldCommonBase::changed() const
{ return changed_during_set_time; }


FieldCommonBase::~FieldCommonBase() {}



/****************************************************************************
 *  Implementation of MultiField
 */




/****************************************************************************
 *  Instances of field templates
 */



INSTANCE_ALL(FieldBase)
INSTANCE_ALL(Field)
INSTANCE_ALL(BCField)
INSTANCE_ALL(MultiField)

