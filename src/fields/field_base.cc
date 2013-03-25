/*
 * field_base.cc
 *
 *  Created on: Dec 4, 2012
 *      Author: jb
 */



#include "fields/field_base_impl.hh"


/****************************************************************************
 *  Implementation of FieldCommon
 */

FieldCommonBase::FieldCommonBase(bool bc)
: n_comp_(0),
  bc_(bc),
  element_selection_(NULL),
  default_( IT::Default::obligatory()),
  mesh_(NULL),
  changed_during_set_time(true), // safe value, we do not miss possible update
  changed_from_last_set_time_(false)
{}

// setters
void FieldCommonBase::set_name(const string & name)
{ name_ =name; }
void FieldCommonBase::set_desc(const string & desc)
{ desc_=desc; }
void FieldCommonBase::set_default(const IT::Default &dflt)
{ default_=dflt;}
void FieldCommonBase::set_n_comp( unsigned int n_comp)
{ n_comp_=n_comp; }
void FieldCommonBase::set_selection( Input::Type::Selection *element_selection)
{ element_selection_=element_selection;}
void FieldCommonBase::set_mesh(Mesh *mesh)
{ mesh_=mesh; }



// getters
const std::string & FieldCommonBase::name() const
{ return name_; }
const std::string & FieldCommonBase::desc() const
{ return desc_; }
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
 *  Instances of field templates
 */



INSTANCE_ALL(FieldBase)
INSTANCE_ALL(Field)
INSTANCE_ALL(BCField)

