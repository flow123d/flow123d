/*
 * field_base.cc
 *
 *  Created on: Dec 4, 2012
 *      Author: jb
 */



#include "fields/field_base_impl.hh"
#include "fields/field_python_impl.hh"
#include "fields/field_constant_impl.hh"
#include "fields/field_formula_impl.hh"
#include "fields/field_interpolated_p0_impl.hh"
#include "fields/field_add_potential_impl.hh"
#include "fields/field_elementwise_impl.hh"




/****************************************************************************
 *  Implementation of FieldCommon
 */

FieldCommonBase::FieldCommonBase(bool bc)
: n_comp_(0),
  bc_(bc),
  element_selection_(NULL),
  default_( IT::Default::obligatory()),
  mesh_(NULL)
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


FieldCommonBase::~FieldCommonBase() {}



/****************************************************************************
 *  Instances of field templates
 */


#define INSTANCE_DIM_DEP_VALUES( field, dim_from, dim_to)                                                               \
template class field<dim_from, FieldValue<dim_to>::VectorFixed >;                       \
template class field<dim_from, FieldValue<dim_to>::TensorFixed >;                       \

//first dimension independent values then dimension dependent
#define INSTANCE_TO_ALL(field, dim_from) \
template class field<dim_from, FieldValue<0>::Enum >;                       \
template class field<dim_from, FieldValue<0>::EnumVector >;                \
template class field<dim_from, FieldValue<0>::Integer >;                       \
template class field<dim_from, FieldValue<0>::Scalar >;                       \
template class field<dim_from, FieldValue<0>::Vector >;                         \
\
INSTANCE_DIM_DEP_VALUES( field, dim_from, 2) \
INSTANCE_DIM_DEP_VALUES( field, dim_from, 3) \

#define INSTANCE_ALL(field) \
INSTANCE_TO_ALL(field, 0) \
INSTANCE_TO_ALL( field, 1) \
INSTANCE_TO_ALL( field, 2) \
INSTANCE_TO_ALL( field, 3)


INSTANCE_ALL(FieldBase)
INSTANCE_ALL(Field)
INSTANCE_ALL(BCField)
INSTANCE_ALL(FieldConstant)
INSTANCE_ALL(FieldPython)
INSTANCE_ALL(FieldFormula)
INSTANCE_ALL(FieldElementwise)
//INSTANCE_ALL(FieldInterpolatedP0)

template class FieldAddPotential<3, FieldValue<0>::Scalar >;
template class FieldAddPotential<2, FieldValue<0>::Scalar >;


