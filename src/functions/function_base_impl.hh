/*
 * function_base_impl.hh
 *
 *  Created on: Oct 1, 2012
 *      Author: jb
 */

#include <functions/functions_all.hh>

#ifndef FUNCTION_BASE_IMPL_HH_
#define FUNCTION_BASE_IMPL_HH_


#include "functions/function_base.hh"
#include "functions/function_interpolated_p0.hh"
#include "functions/function_python.hh"
#include "input/input_type.hh"


template <int dim>
FunctionBase<dim>::FunctionBase(const unsigned int n_components, const double init_time)
: n_components_(n_components), time_(init_time)
{}



template <int dim>
Input::Type::AbstractRecord &FunctionBase<dim>::get_input_type() {
    using namespace Input::Type;
    static AbstractRecord rec("Function", "Abstract record for all time-space functions.");

    if (! rec.is_finished()) {
        rec.finish();
#ifdef HAVE_PYTHON
        FunctionPython<dim>::get_input_type();
#endif
        FunctionInterpolatedP0<dim>::get_input_type();
        rec.no_more_descendants();
    }
    return rec;
}



template <int dim>
FunctionBase<dim> *  FunctionBase<dim>::function_factory(Input::AbstractRecord &rec) {

    if (rec.type() == FunctionInterpolatedP0<dim>::get_input_type()) {
        return  new FunctionInterpolatedP0<dim>(rec);
#ifdef HAVE_PYTHON
    } else if (rec.type() == FunctionPython<dim>::get_input_type()) {
        return  new FunctionPython<dim>(rec);
#endif
    } else {
        xprintf(PrgErr,"TYPE of Function is out of set of descendants. SHOULD NOT HAPPEN.\n");
    }
}


template <int dim>
void FunctionBase<dim>::init_from_input(Input::Record &rec) {
    xprintf(PrgErr, "The function do not support initialization from input.\n");
}



template <int dim>
void FunctionBase<dim>::set_time(double time) {
    time_ = time;
}




template <int dim>
void FunctionBase<dim>::set_element(Element *element) {
    element_ = element;
}

#endif //FUNCTION_BASE_IMPL_HH_
