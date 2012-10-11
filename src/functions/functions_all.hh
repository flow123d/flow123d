/*
 * functions_all.hh
 *
 *  Created on: Oct 10, 2012
 *      Author: jb
 */

#ifndef FUNCTIONS_ALL_HH_
#define FUNCTIONS_ALL_HH_

/**
 * This is temporary solution for the problem that comes with Input::Type::AbstractRecord.
 * the static function that construct the AbstractRecord for functions in abstract class template  FunctionBase<dim>
 * has to call get_input_type static methods of its descendants. To this end we need first definition of all classes and then
 * definition of its methods.
 *
 * The true solution is in better way how to construct Input types. There should be way how to construct an input object
 * on one line (like in boost::program_options), then we can replace static methods with static data members (that are initialized at
 * very beginning of the program) and than it is enough for descendants of an AbstractRecord to "register" to
 * its parent.
 *
 * To guarantee correct order of static initializations, one can use :
 *
 * http://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/Nifty_Counter
 */

#include <functions/function_base.hh>
#include <functions/function_interpolated_p0.hh>
#include <functions/function_python.hh>
#include <functions/function_elementwise.hh>

#include <functions/function_base_impl.hh>
#include <functions/function_interpolated_p0_impl.hh>
#include <functions/function_python.hh>
#include <functions/function_elementwise_impl.hh>

#endif /* FUNCTIONS_ALL_HH_ */
