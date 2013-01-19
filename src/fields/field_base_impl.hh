/*
 * function_base_impl.hh
 *
 *  Created on: Oct 1, 2012
 *      Author: jb
 */


#ifndef FIELD_BASE_IMPL_HH_
#define FIELD_BASE_IMPL_HH_

#include <string>
using namespace std;

#include <boost/type_traits.hpp>
#include <boost/format.hpp>


#include "fields/field_base.hh"
#include "fields/field_interpolated_p0.hh"
#include "fields/field_python.hh"
#include "fields/field_constant.hh"
#include "fields/field_formula.hh"

#include "fields/field_values.hh"

#include "input/input_type.hh"


namespace it = Input::Type;







// ************************************************************************************************** implementation of FieldBase<...>

template <int spacedim, class Value>
it::AbstractRecord FieldBase<spacedim, Value>::input_type
    = it::AbstractRecord("Field"+template_name(), "Abstract record for all time-space functions.")
          .allow_auto_conversion("FieldConstant");


template <int spacedim, class Value>
FieldBase<spacedim, Value>::FieldBase(const double init_time, unsigned int n_comp)
: time_(init_time), value_(r_value_)
{
    value_.set_n_comp(n_comp);
}




template <int spacedim, class Value>
string FieldBase<spacedim, Value>::template_name() {
    return boost::str(boost::format("R%i -> %s") % spacedim % Value::type_name() );
}


template <int spacedim, class Value>
Input::Type::AbstractRecord FieldBase<spacedim, Value>::get_input_type(typename Value::ElementInputType *element_input_type) {
    // jak sem dostat primo potomky
    it::AbstractRecord type= it::AbstractRecord("Field", "Abstract record for all time-space functions.");
    type.allow_auto_conversion("FieldConstant");

    FieldConstant<spacedim,Value>::get_input_type(type, element_input_type);
    FieldFormula<spacedim,Value>::get_input_type(type, element_input_type);
#ifdef HAVE_PYTHON
    FieldPython<spacedim,Value>::get_input_type(type, element_input_type);
#endif
    //FieldInterpolatedP0<spacedim,Value>::get_input_type(type, element_input_type);

    return type;
}



template <int spacedim, class Value>
FieldBase<spacedim, Value> *  FieldBase<spacedim, Value>::function_factory(
        Input::AbstractRecord rec, double init_time, unsigned int n_comp )
{
    FieldBase<spacedim, Value> *func;

    if (rec.type() == FieldInterpolatedP0<spacedim,Value>::input_type ) {
//        func= new FieldInterpolatedP0<spacedim,Value>(init_time, n_comp);
#ifdef HAVE_PYTHON
    } else if (rec.type() == FieldPython<spacedim,Value>::input_type ) {
        func= new FieldPython<spacedim, Value>(init_time, n_comp);
#endif
    } else if (rec.type() == FieldConstant<spacedim, Value>::input_type ) {
        func=new FieldConstant<spacedim,Value>(init_time, n_comp);
    } else if (rec.type() == FieldFormula<spacedim,Value>::input_type ) {
        func=new FieldFormula<spacedim,Value>(init_time, n_comp);
    } else {
        xprintf(PrgErr,"TYPE of Field is out of set of descendants. SHOULD NOT HAPPEN.\n");
    }
    func->init_from_input(rec);
    return func;
}



template <int spacedim, class Value>
void FieldBase<spacedim, Value>::init_from_input(Input::Record rec) {
    xprintf(PrgErr, "The function do not support initialization from input.\n");
}



template <int spacedim, class Value>
void FieldBase<spacedim, Value>::set_time(double time) {
    time_ = time;
}


/*
template <int spacedim, class Value>
typename Value::return_type &FieldBase<spacedim, Value>::value(const Point<spacedim> &p, ElementAccessor<spacedim> &elm) {
    this->value(p,elm, this->r_value_);
    return this->r_value_;
}
*/

/*
template <int spacedim, class Value>
void FieldBase<spacedim, Value>::value_list(const std::vector< Point<spacedim> >  &point_list, ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type>  &value_list, std::vector<FieldResult> &result_list)
{
    ASSERT_SIZES( point_list.size(), value_list.size() );
    for(unsigned int i=0; i< point_list.size(); i++)
        value_list[i] = value(point_list[i], elm);
}
*/



// ************************************************************************************************** implementation of Field<...>


template<int spacedim, class Value>
Field<spacedim,Value>::Field()
: FieldCommonBase(false)
{
    this->enum_valued_ = boost::is_same<typename Value::element_type, FieldEnum>::value;
}



template<int spacedim, class Value>
typename Field<spacedim,Value>::FieldBaseType * Field<spacedim,Value>::operator() (Region reg) {
    ASSERT_LESS(reg.idx(), this->region_fields.size());
    return this->region_fields[reg.idx()];
}



template<int spacedim, class Value>
void Field<spacedim, Value>::init_from_input(Region reg, Input::AbstractRecord rec) {
    // initialize table if it is empty, we assume that the RegionDB is closed at this moment
    if (region_fields.size() == 0)
        region_fields.resize(Region::db().size(), NULL);

    if (region_fields[reg.idx()] != NULL) {
        xprintf(
                Warn, "Overwriting existing value on region ID=%d. In initialization of the Field: '%s'.\n", reg.id(), this->name().c_str());
        delete region_fields[reg.idx()];
    }
    region_fields[reg.idx()] = FieldBaseType::function_factory(rec, 0.0, this->n_comp_);
}


template<int spacedim, class Value>

void Field<spacedim, Value>::set_time(double time) {
    for(unsigned int i=0; i < region_fields.size(); i++) region_fields[i]->set_time(time);
}



template<int spacedim, class Value>
inline typename Value::return_type & Field<spacedim,Value>::value(const Point<spacedim> &p, ElementAccessor<spacedim> &elm)  {
    ASSERT_LESS(elm.region().idx(), region_fields.size() );
    ASSERT( region_fields[elm.region().idx()] , "Null field ptr on region %d, field: %s\n", elm.region().idx(), this->name_.c_str());
    return region_fields[elm.region().idx()]->value(p,elm);
}



template<int spacedim, class Value>
inline FieldResult Field<spacedim,Value>::value_list(const std::vector< Point<spacedim> >  &point_list, ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type>  &value_list)
{
    region_fields[elm.region().idx()]->value_list(point_list,elm, value_list);
}



#endif //FUNCTION_BASE_IMPL_HH_
