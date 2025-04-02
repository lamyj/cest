#include "bm2.h"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include "misc.h"

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

void bm2(pybind11::module & m)
{
    using namespace pybind11::literals;
    
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
