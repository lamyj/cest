#ifndef _dd19598e_515e_45b2_a47a_8068ab5daef9
#define _dd19598e_515e_45b2_a47a_8068ab5daef9

#include <Eigen/Core>

#include <pybind11/pybind11.h>

#include "misc.h"

Eigen::Matrix10d
bm(
    Species const & species_a, Species const & species_b, Species const & species_c,
    double Cb, double Cc,
    double w0, double delta_w_rf, double w1, double duration);

Eigen::Vector9d bm(
    Species const & species_a, Species const & species_b, Species const & species_c,
    double Cb, double Cc,
    double w0, double delta_w_rf, double w1, double duration,
    Eigen::Vector9d const & M);

Eigen::Vector10d
bm(
    Species const & species_a, Species const & species_b, Species const & species_c,
    double Cb, double Cc,
    double w0, double delta_w_rf, Eigen::VectorXd const & w1,
    double step, Eigen::Vector10d const & M0);

Eigen::Vector9d
bm(
    Species const & species_a, Species const & species_b, Species const & species_c,
    double Cb, double Cc,
    double w0, double delta_w_rf, Eigen::VectorXd const & w1,
    double step, Eigen::Vector9d const & M0);

void bm3_partial(pybind11::module & m);

#endif // _dd19598e_515e_45b2_a47a_8068ab5daef9
