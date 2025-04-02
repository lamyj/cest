#ifndef _0b16c428_bf5f_41fe_bd25_64831665ce8a
#define _0b16c428_bf5f_41fe_bd25_64831665ce8a

#include <Eigen/Core>

#include <pybind11/pybind11.h>

#include "misc.h"

/**
 * @brief Two-pools Bloch-McConnell transition matrix
 * @param species_a
 * @param species_b
 * @param Cb Transition rate from B to A (Hz)
 * @param w0 Larmor frequency (rad/s)
 * @param delta_w_rf Frequency offset of the saturation RF pulse (ppm)
 * @param w1 Frequency of the B1 field of the saturation RF pulse (rad/s)
 * @param duration Duration of the simulation in s
 * @return 7x7 matrix in projective coordinates, should be applied to a vector
 * [Mxa, Mya, Mza, Mxb, Myb, Mzb, 1]
 */
Eigen::Matrix7d bm(
    Species const & species_a, Species const & species_b, double Cb,
    double w0, double delta_w_rf, double w1, double duration);

/**
 * @brief Two-pools Bloch-McConnell simulation
 * @param species_a
 * @param species_b
 * @param Cb Transition rate from B to A (Hz)
 * @param w0 Larmor frequency (rad/s)
 * @param delta_w_rf Frequency offset of the saturation RF pulse (ppm)
 * @param w1 Frequency of the B1 field of the saturation RF pulse (rad/s)
 * @param duration Duration of the simulation in s
 * @param M magnetization as [Mxa, Mya, Mza, Mxb, Myb, Mzb]
 * @return Magnetization after evolution, as [Mxa, Mya, Mza, Mxb, Myb, Mzb]
 * 
 */
Eigen::Vector6d bm(
    Species const & species_a, Species const & species_b, double Cb,
    double w0, double delta_w_rf, double w1, double duration,
    Eigen::Vector6d const & M);

/**
 * @brief Two-pools Bloch-McConnel simulation of a shaped pulse
 * @param species_a
 * @param species_b
 * @param Cb Transition rate from B to A (Hz)
 * @param w0 Larmor frequency (rad/s)
 * @param delta_w_rf Frequency offset of the saturation RF pulse (ppm)
 * @param w1 Frequency of the B1 field of the saturation RF pulse (rad/s)
 * @param step Interval between to w1 values of the simulation in s
 * @param M0 Initial magnetization as [Mxa, Mya, Mza, Mxb, Myb, Mzb, 1]
 * @return Magnetization after evolution, as [Mxa, Mya, Mza, Mxb, Myb, Mzb, 1]
 */
Eigen::Vector7d bm(
    Species const & species_a, Species const & species_b, double Cb,
    double w0, double delta_w_rf, Eigen::VectorXd const & w1,
    double step, Eigen::Vector7d const & M0);

/**
 * @brief Two-pools Bloch-McConnel simulation of a shaped pulse
 * @param species_a
 * @param species_b
 * @param Cb Transition rate from B to A (Hz)
 * @param w0 Larmor frequency (rad/s)
 * @param delta_w_rf Frequency offset of the saturation RF pulse (ppm)
 * @param w1 Frequency of the B1 field of the saturation RF pulse (rad/s)
 * @param step Interval between to w1 values of the simulation in s
 * @param M0 Initial magnetization as [Mxa, Mya, Mza, Mxb, Myb, Mzb]
 * @return Magnetization after evolution, as [Mxa, Mya, Mza, Mxb, Myb, Mzb]
 */
Eigen::Vector6d
bm(
    Species const & species_a, Species const & species_b, double Cb,
    double w0, double delta_w_rf, Eigen::VectorXd const & w1,
    double step, Eigen::Vector6d const & M0);

void bm2(pybind11::module & m);

#endif // _0b16c428_bf5f_41fe_bd25_64831665ce8a
