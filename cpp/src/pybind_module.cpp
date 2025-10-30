#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(_sirius, m) {
    m.doc() = "Tiny module (diagnostics)";
    m.def("ping", [](){ return "pong"; });
}
