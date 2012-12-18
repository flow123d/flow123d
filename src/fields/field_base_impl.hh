/*
 * function_base_impl.hh
 *
 *  Created on: Oct 1, 2012
 *      Author: jb
 */


#include <fields/field_all.hh>

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
#include "fields/field_values.hh"

#include "input/input_type.hh"


namespace it = Input::Type;

template <int spacedim, class Value>
it::AbstractRecord FieldBase<spacedim, Value>::input_type
    = it::AbstractRecord("Field", "Abstract record for all time-space functions.");


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
FieldBase<spacedim, Value> *  FieldBase<spacedim, Value>::function_factory(
        Input::AbstractRecord rec, double init_time, unsigned int n_comp )
{
    FieldBase<spacedim, Value> *func;

    if (rec.type() == FieldInterpolatedP0<spacedim,Value>::get_input_type()) {
//        func= new FieldInterpolatedP0<spacedim,Value>(init_time, n_comp);
#ifdef HAVE_PYTHON
    } else if (rec.type() == FieldPython<spacedim,Value>::get_input_type()) {
        func= new FieldPython<spacedim, Value>(init_time, n_comp);
#endif
    } else if (rec.type() == FieldConstant<spacedim, Value>::get_input_type()) {
        func=new FieldConstant<spacedim,Value>(init_time, n_comp);
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



template <int spacedim, class Value>
typename Value::return_type &FieldBase<spacedim, Value>::value(const Point<spacedim> &p, ElementAccessor<spacedim> &elm) {
    this->value(p,elm, this->r_value_);
    return this->r_value_;
}



template <int spacedim, class Value>
void FieldBase<spacedim, Value>::value_list(const std::vector< Point<spacedim> >  &point_list, ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type>  &value_list, std::vector<FieldResult> &result_list)
{
    ASSERT_SIZES( point_list.size(), value_list.size() );
    for(unsigned int i=0; i< point_list.size(); i++)
        result_list[i] = value(point_list[i], elm, value_list[i]);
}

#endif //FUNCTION_BASE_IMPL_HH_
