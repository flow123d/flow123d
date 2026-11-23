/**
 * @file python_wrapper.cc
 * @brief Python entry point for running Flow123d.
 */

#include <pybind11/pybind11.h>

#include "coupling/application.hh"

PYBIND11_MODULE(flow123d_python_api, m) {
    m.doc() = "Python interface to Flow123d";

    m.def(
        "run",
        [](const pybind11::dict &input, const pybind11::dict &env_args) {
            Application app;
            app.init_python(env_args);

            Input::Record record = app.read_input(input);
            app.run(record);
        },
        pybind11::arg("input"),
		pybind11::arg("env_args") = pybind11::dict()
    );
}
