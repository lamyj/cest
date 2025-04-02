#include "misc.h"

#include <Eigen/Core>

#include <pybind11/pybind11.h>

void misc(pybind11::module m)
{
    using namespace pybind11::literals;
    
    pybind11::class_<Species>(m, "Species", "Chemical species")
        .def_readwrite("T1", &Species::T1, "T1 relaxation times in s")
        .def_readwrite("T2", &Species::T2, "T2 relaxation times in s")
        .def_readwrite("delta_w", &Species::delta_w, "Frequency offset in ppm")
        .def_readwrite("M0", &Species::M0, "Equilibrium magnetization (unitless)")
        .def(
            pybind11::init(
                [](double T1, double T2, double delta_w, double M0){
                    return Species{T1, T2, delta_w, M0};
                }),
            "T1"_a, "T2"_a, "delta_w"_a, "M0"_a)
        .def(pybind11::pickle(
            [](Species const & species) {
                pybind11::dict state;
                state["T1"] = species.T1;
                state["T2"] = species.T2;
                state["delta_w"] = species.delta_w;
                state["M0"] = species.M0;
                return state;
            },
            [](pybind11::dict state) {
                Species species;
                species.T1 = state["T1"].cast<double>();
                species.T2 = state["T2"].cast<double>();
                species.delta_w = state["delta_w"].cast<double>();
                species.M0 = state["M0"].cast<double>();
                return species;
            }));
}
