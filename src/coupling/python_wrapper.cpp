/**
 * @file python_wrapper.cpp
 * @brief Python entry point for running Flow123d.
 */

#include <pybind11/pybind11.h>

#include "coupling/application.hh"

PYBIND11_MODULE(flow123d_main_api, m) {
    m.doc() = "Python interface to Flow123d";

    m.def(
        "run",
        [](const pybind11::dict &input) {
            Application app;
            app.init_python();

            Input::Record record = app.read_input(input);
            app.run(record);
        },
        pybind11::arg("input")
    );
}
