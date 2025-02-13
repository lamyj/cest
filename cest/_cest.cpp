#include <cmath>

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace Eigen
{
using Matrix6d = Eigen::Matrix<double, 6, 6>;
using Vector6d = Eigen::Matrix<double, 6, 1>;
using RowVector6d = Eigen::Matrix<double, 1, 6>;

using Matrix7d = Eigen::Matrix<double, 7, 7>;
using Vector7d = Eigen::Matrix<double, 7, 1>;
using RowVector7d = Eigen::Matrix<double, 1, 7>;
}

/// @brief Chemical species
class Species
{
public:
    /// @brief T1 relaxation times in s
    double T1;
    /// @brief T2 relaxation times in s
    double T2;
    /// @brief Frequency offset in ppm
    double delta_w;
    /// @brief Equilibrium magnetization (unitless)
    double M0;
};

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
    double w0, double delta_w_rf, double w1, double duration)
{
    auto const M0a = species_a.M0;
    auto const R1a = 1/species_a.T1;
    auto const R2a = 1/species_a.T2;
    
    auto const M0b = species_b.M0;
    auto const R1b = 1/species_b.T1;
    auto const R2b = 1/species_b.T2;
    
    auto const Ca = M0b/M0a * Cb;
    
    auto const k1a = R1a+Ca;
    auto const k2a = R2a+Ca;
    
    auto const k1b = R1b+Cb;
    auto const k2b = R2b+Cb;
    
    auto const wa = w0*(1 + species_a.delta_w*1e-6);
    auto const wb = w0*(1 + species_b.delta_w*1e-6);
    
    auto const w = w0*(1 + delta_w_rf*1e-6);
    
    Eigen::Matrix7d dM;
    dM << 
    //   Mxa      Mya   Mza   Mxa      Myb   Mzb            M0
        -k2a, -(wa-w),    0,   Cb,       0,    0,            0,
        wa-w,    -k2a,  -w1,    0,      Cb,    0,            0,
           0,      w1, -k1a,    0,       0,   Cb,      M0a*R1a,
          Ca,       0,    0, -k2b, -(wb-w),    0,            0,
           0,      Ca,    0, wb-w,    -k2b,  -w1,            0,
           0,       0,   Ca,    0,      w1, -k1b,      M0b*R1b,
           0,       0,    0,    0,       0,    0,            0;
    return (dM*duration).exp();
}

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
    Eigen::Vector6d const & M)
{
    auto const M0a = species_a.M0;
    auto const R1a = 1/species_a.T1;
    auto const R2a = 1/species_a.T2;
    
    auto const M0b = species_b.M0;
    auto const R1b = 1/species_b.T1;
    auto const R2b = 1/species_b.T2;
    
    auto const Ca = M0b/M0a * Cb;
    
    auto const k1a = R1a+Ca;
    auto const k2a = R2a+Ca;
    
    auto const k1b = R1b+Cb;
    auto const k2b = R2b+Cb;
    
    auto const wa = w0*(1 + species_a.delta_w*1e-6);
    auto const wb = w0*(1 + species_b.delta_w*1e-6);
    
    auto const w = w0*(1 + delta_w_rf*1e-6);
    
    Eigen::Matrix6d A;
    A << 
    //   Mxa      Mya   Mza   Mxa      Myb   Mzb
        -k2a, -(wa-w),    0,   Cb,       0,    0,
        wa-w,    -k2a,  -w1,    0,      Cb,    0,
           0,      w1, -k1a,    0,       0,   Cb,
          Ca,       0,    0, -k2b, -(wb-w),    0,
           0,      Ca,    0, wb-w,    -k2b,  -w1,
           0,       0,   Ca,    0,      w1, -k1b;
    
    Eigen::Vector6d const b{0, 0, M0a*R1a, 0, 0, M0b*R1b};
    
    Eigen::Vector6d const AinvB = A.partialPivLu().solve(b);
    
    return (A*duration).exp() * (M+AinvB) - AinvB;
}

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
Eigen::Vector7d
bm(
    Species const & species_a, Species const & species_b, double Cb,
    double w0, double delta_w_rf, Eigen::VectorXd const & w1,
    double step, Eigen::Vector7d const & M0)
{
    Eigen::Vector7d M = M0;
    for(auto && w1_: w1)
    {
        M = bm(species_a, species_b, Cb, w0, delta_w_rf, w1_, step) * M;
    }
    return M;
}

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
    double step, Eigen::Vector6d const & M0)
{
    Eigen::Vector6d M = M0;
    for(auto && w1_: w1)
    {
        M = bm(species_a, species_b, Cb, w0, delta_w_rf, w1_, step, M);
    }
    return M;
}

PYBIND11_MODULE(_cest, m)
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
            "T1"_a, "T2"_a, "delta_w"_a, "M0"_a);
    
    m.def(
        "bm",
        pybind11::overload_cast<
            Species const &, Species const &,
            double, double, double, double, double>(&bm),
        "Two-pools Bloch-McConnell transition matrix between two pools A and B, "
        "using projective coordinates",
        "species_a"_a, "species_b"_a,
        "Cb"_a, "w0"_a, "delta_w_rf"_a, "w1"_a, "step"_a);
    
    m.def(
        "bm",
        pybind11::overload_cast<
            Species const &, Species const &,
            double, double, double, double, double, Eigen::Vector6d const &>(&bm),
        "Two-pools Bloch-McConnell simulation",
        "species_a"_a, "species_b"_a, "Cb"_a, "w0"_a, "delta_w_rf"_a, "w1"_a,
        "step"_a, "M0"_a);
    
    m.def(
        "bm",
        pybind11::overload_cast<
            Species const &, Species const &,
            double, double, double, Eigen::VectorXd const &, double,
            Eigen::Vector7d const &>(&bm),
        "Two-pools Bloch-McConnel simulation of a shaped pulse",
        "species_a"_a, "species_b"_a, "Cb"_a, "w0"_a, "delta_w_rf"_a, "w1"_a,
        "step"_a, "M0"_a);
    
    m.def(
        "bm",
        pybind11::overload_cast<
            Species const &, Species const &,
            double, double, double, Eigen::VectorXd const &, double,
            Eigen::Vector6d const &>(&bm),
        "Two-pools Bloch-McConnel simulation of a shaped pulse",
        "species_a"_a, "species_b"_a, "Cb"_a, "w0"_a, "delta_w_rf"_a, "w1"_a,
        "step"_a, "M0"_a);
}
