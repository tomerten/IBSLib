#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
using namespace pybind11::literals;

#include "NumericFunctions/NumericFunctions.hpp"
#include "Twiss/Twiss.hpp"

PYBIND11_MODULE(ibslib_pb, m) {
  m.doc() = "Twiss reader cpp";
  m.def("GetTwissHeader", &GetTwissHeader,
        "Get the twiss header as a dictionary.");
  m.def("GetTwissTable", &GetTwissTableAsMap, "Get the twiss data table.");
  m.attr("clight") = py::float_(clight);
  m.attr("hbarGeV") = py::float_(hbar);
  m.attr("electron_mass") = py::float_(emass);
  m.attr("proton_mass") = py::float_(pmass);
  m.attr("neutron_mass") = py::float_(nmass);
  m.attr("mu_mass") = py::float_(mumass);
  m.attr("atomic_mass_unit") = py::float_(atomicmassunit);
  m.attr("pi") = py::float_(pi);
  m.attr("electric_charge") = py::float_(ec);
  m.attr("euler") = py::float_(euler);
  m.attr("electron_radius") = py::float_(erad);
  m.attr("proton_radius") = py::float_(prad);
  m.def("sige_from_sigs", &sigefromsigs, "Energy spread from bunch length.");
  m.def("sigs_from_sige", &sigsfromsige, "Bunch length from energy spread.");
  /*
  m.def("eta", &eta, "Slip factor");
  m.def("fmohl", &fmohl, "Fmohl function");
  m.def("particle_radius", &ParticleRadius,
        "Particle radius from charge and atomic mass.");
  m.def("beta_relativistic_from_gamma", &BetaRelativisticFromGamma,
        "Beta relativistic.");
   m.def("rds", &rds, "Nagaitsev paper rds function");*/
}