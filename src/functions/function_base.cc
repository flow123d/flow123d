/*
 * function_base.cc
 *
 *  Created on: Oct 1, 2012
 *      Author: jb
 */

#include "functions/function_base.hh"
#include "input/input_type.hh"

Input::Type::AbstractRecord &FunctionBase::get_input_type() {
    using namespace Input::Type;
    static AbstractRecord rec("Function", "Abstract record for all time-space functions.");

    if (! rec.finished()) {
        rec.finish();

        FunctionPython::get_input_type();
        FunctionElementwise::get_input_type();
        rec.no_more_descendants();
    }
    return rec;
}




