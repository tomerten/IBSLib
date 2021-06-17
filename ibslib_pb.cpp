#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
using namespace pybind11::literals;

#include "Twiss/Twiss.hpp"

PYBIND11_MODULE(ibslib_pb, m) {
  m.doc() = "Twiss reader cpp";
  m.def("GetTwissHeader", &GetTwissHeader, "Get the twiss header as a map.");
}