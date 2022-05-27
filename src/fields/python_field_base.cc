/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 *
 * @file    python_field_base.cc
 * @brief
 */

#include <include/pybind11/pybind11.h>
#include <include/pybind11/embed.h>    // everything needed for embedding
#include <include/pybind11/stl.h>      // type conversion
#include "fields/python_field_base.hh"

//int add(int i, int j) {
//    return i + j;
//}

PYBIND11_MODULE(flowpy, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring
//    m.def("add", &add, "A function that adds two numbers");

    py::class_<PythonFieldBase>(m, "PythonFieldBase")
        .def(py::init<>())
        //.def(py::init<std::string, std::vector<ssize_t>, std::vector<double>>())
        .def("set_result", &PythonFieldBase::set_result)
		.def("add_to_dict", &PythonFieldBase::add_to_dict)
        .def("set_result_data", &PythonFieldBase::set_result_data)
		.def("add_to_dict_data", &PythonFieldBase::add_to_dict_data)
		.def("print_fields", &PythonFieldBase::print_fields)
		.def("print_result", &PythonFieldBase::print_result)
        .def_property("t", &PythonFieldBase::get_time, &PythonFieldBase::set_time)
	    .def_property("result", &PythonFieldBase::get_field_result, &PythonFieldBase::set_field_result)
	    .def_property("f_dict", &PythonFieldBase::get_fields_dict, &PythonFieldBase::set_fields_dict);

}

//void fill_vec(std::vector<double> &vec) {
//	for (uint i=0; i<vec.size(); ++i) vec[i] = 1.0 + i;
//}
//
//int main() {
//	py::scoped_interpreter guard{}; // start the interpreter and keep it alive
//
//	// test of simple function in C++
//	int a = 1;
//	int b = 2;
//	int c = add(a, b);
//	std::cout << "Add function: " << a << " + " << b << " = " << c << std::endl;
//
//	// test of PythonFieldBase object in C++
//	std::vector<double> csection_vec(16);
//	fill_vec(csection_vec);
//	std::vector<ssize_t> csection_shape = {1,1};
//	std::vector<double> velocity_vec(48);
//	fill_vec(velocity_vec);
//	std::vector<ssize_t> velocity_shape = {1,3};
//	std::vector<double> result_vec(48);
//	fill_vec(result_vec);
//	std::vector<ssize_t> result_shape = {1,3};
//	std::vector<FieldCacheProxy> field_data;
//	field_data.emplace_back("csection", csection_shape, csection_vec);
//	field_data.emplace_back("velocity", velocity_shape, velocity_vec);
//	FieldCacheProxy result_data("result", result_shape, result_vec);
//	PythonFieldBase field(field_data, result_data);
//	field.print_fields();
//
//	// test of call of simple function in Python
//	// source: https://stackoverflow.com/questions/42521830/call-a-python-function-from-c-using-pybind11
//	auto math = py::module::import("math");
//	double root_two = math.attr("sqrt")(2.0).cast<double>();
//	std::cout << "The square root of 2 is: " << root_two << "\n";
//
//	// test of call of 'add' function in Python
//	// source: same as previous
//	auto flowpy = py::module::import("flowpy");
//    py::function add_func =
//        py::reinterpret_borrow<py::function>(   // cast from 'object' to 'function - use `borrow` (copy) or `steal` (move)
//            py::module::import("flowpy").attr("add")  // import method "min_rosen" from python "module"
//        );
//	int add_result = add_func(2, 3).cast<int>();
//	std::cout << "Result of add(2, 3) is: " << add_result << "\n";
//	std::cout << "-----------------------------------------------------\n";
//
//	// test of using function - use object 'PythonFieldBase field;'
//	std::cout << "The most important test simulated behavior of Flow123D new assembly.\n";
//	py::function fn_init =
//	    py::reinterpret_borrow<py::function>(   // cast from 'object' to 'function - use `borrow` (copy) or `steal` (move)
//	        py::module::import("example_func").attr("fn_init")  // import method "min_rosen" from python "module"
//	    );
//	py::function fn_reinit =
//	    py::reinterpret_borrow<py::function>(   // cast from 'object' to 'function - use `borrow` (copy) or `steal` (move)
//	        py::module::import("example_func").attr("fn_reinit")  // import method "min_rosen" from python "module"
//	    );
//	py::function fn_eval =
//	    py::reinterpret_borrow<py::function>(   // cast from 'object' to 'function - use `borrow` (copy) or `steal` (move)
//	        py::module::import("example_func").attr("fn_eval")  // import method "min_rosen" from python "module"
//	    );
//	py::object user_context = fn_init(field);
//	py::object f = fn_reinit(field, user_context);
//	fn_eval(field, user_context);
//	std::cout << "After 'fn_eval' function:" << std::endl;
//	field.print_result();
//
//	//https://pybind11.readthedocs.io/en/stable/advanced/pycpp/object.html#calling-python-methods
//
//	return 0;
//}
