/*
 * function_base_impl.hh
 *
 *  Created on: Oct 1, 2012
 *      Author: jb
 */


#ifndef field_algo_base_IMPL_HH_
#define field_algo_base_IMPL_HH_

#include <string>
#include <limits>
#include <memory>
using namespace std;

#include "fields/field_algo_base.hh"
#include "fields/field_interpolated_p0.hh"
#include "fields/field_python.hh"
#include "fields/field_constant.hh"
#include "fields/field_formula.hh"
#include "fields/field_elementwise.hh"

#include "fields/field_values.hh"

#include "input/factory.hh"


namespace it = Input::Type;


/******************************************************************************************
 * Implementation of FieldBase<...>
 */

template <int spacedim, class Value>
FieldAlgorithmBase<spacedim, Value>::FieldAlgorithmBase(unsigned int n_comp)
: time_( -numeric_limits<double>::infinity() ),
  value_(r_value_)
{
    value_.set_n_comp(n_comp);
}



template <int spacedim, class Value>
string FieldAlgorithmBase<spacedim, Value>::template_name() {
    return boost::str(boost::format("R%i -> %s") % spacedim % Value::type_name() );
}



template <int spacedim, class Value>
it::AbstractRecord FieldAlgorithmBase<spacedim, Value>::input_type
    = it::AbstractRecord("Field:"+template_name(), "Abstract record for all time-space functions.")
          .allow_auto_conversion("FieldConstant");



template <int spacedim, class Value>
Input::Type::AbstractRecord FieldAlgorithmBase<spacedim, Value>::get_input_type(const typename Value::ElementInputType *element_input_type) {
    it::AbstractRecord type= it::AbstractRecord("Field:"+template_name(), "Abstract record for all time-space functions.");
    type.allow_auto_conversion("FieldConstant");

    FieldConstant<spacedim,Value>::get_input_type(type, element_input_type);
    FieldFormula<spacedim,Value>::get_input_type(type, element_input_type);
#ifdef HAVE_PYTHON
    FieldPython<spacedim,Value>::get_input_type(type, element_input_type);
#endif
    FieldInterpolatedP0<spacedim,Value>::get_input_type(type, element_input_type);
    FieldElementwise<spacedim,Value>::get_input_type(type, element_input_type);

    return type;
}



template <int spacedim, class Value>
shared_ptr< FieldAlgorithmBase<spacedim, Value> >
FieldAlgorithmBase<spacedim, Value>::function_factory(const Input::AbstractRecord &rec, unsigned int n_comp )
{
    shared_ptr< FieldAlgorithmBase<spacedim, Value> > func;
    func = rec.factory< FieldAlgorithmBase<spacedim, Value> >(n_comp);
    func->init_from_input(rec);
    return func;
}



template <int spacedim, class Value>
void FieldAlgorithmBase<spacedim, Value>::init_from_input(const Input::Record &rec) {
    xprintf(PrgErr, "The field '%s' do not support initialization from input.\n",
            typeid(this).name());
}



template <int spacedim, class Value>
bool FieldAlgorithmBase<spacedim, Value>::set_time(double time) {
    time_ = time;
    return false; // no change
}



template <int spacedim, class Value>
void FieldAlgorithmBase<spacedim, Value>::set_mesh(const Mesh *mesh,  bool boundary_domain) {
}



template<int spacedim, class Value>
unsigned int FieldAlgorithmBase<spacedim, Value>::n_comp() const {
    return (Value::NRows_ ? 0 : value_.n_rows());
}



template<int spacedim, class Value>
void FieldAlgorithmBase<spacedim, Value>::value_list(
        const std::vector< Point >  &point_list,
        const ElementAccessor<spacedim> &elm,
        std::vector<typename Value::return_type>  &value_list)
{
    ASSERT_EQUAL( point_list.size(), value_list.size() );
    for(unsigned int i=0; i< point_list.size(); i++) {
        ASSERT( Value(value_list[i]).n_rows()==this->value_.n_rows(),
                "value_list[%d] has wrong number of rows: %d; should match number of components: %d\n",
                i, Value(value_list[i]).n_rows(),this->value_.n_rows());
        value_list[i]=this->value(point_list[i], elm);
    }

}



/****************************************************************************
 *  Macros for explicit instantiation of particular field class template.
 */


// Instantiation of fields with values dependent of the dimension of range space
#define INSTANCE_DIM_DEP_VALUES( field, dim_from, dim_to)                                                               \
template class field<dim_from, FieldValue<dim_to>::VectorFixed >;                       \
template class field<dim_from, FieldValue<dim_to>::TensorFixed >;                       \

// Instantiation of fields with domain in the ambient space of dimension @p dim_from
#define INSTANCE_TO_ALL(field, dim_from) \
template class field<dim_from, FieldValue<0>::Enum >;                       \
template class field<dim_from, FieldValue<0>::EnumVector >;                \
template class field<dim_from, FieldValue<0>::Integer >;                       \
template class field<dim_from, FieldValue<0>::Scalar >;                       \
template class field<dim_from, FieldValue<0>::Vector >;                         \
\
INSTANCE_DIM_DEP_VALUES( field, dim_from, 2) \
INSTANCE_DIM_DEP_VALUES( field, dim_from, 3) \

// All instances of one field class template @p field.
// currently we need only fields on 3D ambient space (and 2D for some tests)
// so this is to save compilation time and avoid memory problems on the test server
#define INSTANCE_ALL(field) \
INSTANCE_TO_ALL( field, 3) \
INSTANCE_TO_ALL( field, 2)
// currently we use only 3D ambient space




#endif //FUNCTION_BASE_IMPL_HH_
