#include "bm3_partial.h"

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include "misc.h"

Eigen::Matrix10d
bm(
    Species const & species_a, Species const & species_b, Species const & species_c,
    double Cb, double Cc,
    double w0, double delta_w_rf, double w1, double duration)
{
    auto const M0a = species_a.M0;
    auto const R1a = 1/species_a.T1;
    auto const R2a = 1/species_a.T2;
    
    auto const M0b = species_b.M0;
    auto const R1b = 1/species_b.T1;
    auto const R2b = 1/species_b.T2;
    
    auto const M0c = species_c.M0;
    auto const R1c = 1/species_c.T1;
    auto const R2c = 1/species_c.T2;
    
    auto const Cab = M0b/M0a * Cb;
    auto const Cac = M0c/M0a * Cc;
    auto const Ca = Cab + Cac;
    
    auto const k1a = R1a+Ca;
    auto const k2a = R2a+Ca;
    
    auto const k1b = R1b+Cb;
    auto const k2b = R2b+Cb;
    
    auto const k1c = R1c+Cc;
    auto const k2c = R2c+Cc;
    
    auto const wa = w0*(1 + species_a.delta_w*1e-6);
    auto const wb = w0*(1 + species_b.delta_w*1e-6);
    auto const wc = w0*(1 + species_b.delta_w*1e-6);
    
    auto const w = w0*(1 + delta_w_rf*1e-6);
    
    Eigen::Matrix10d dM;
    dM << 
    //   Mxa      Mya   Mza   Mxb      Myb   Mzb   Mxc      Myc   Mzc            M0
        -k2a, -(wa-w),    0,   Cb,       0,    0,   Cc,       0,    0,            0,
        wa-w,    -k2a,  -w1,    0,      Cb,    0,    0,      Cc,    0,            0,
           0,      w1, -k1a,    0,       0,   Cb,    0,       0,   Cc,      M0a*R1a,
         Cab,       0,    0, -k2b, -(wb-w),    0,    0,       0,    0,            0,
           0,     Cab,    0, wb-w,    -k2b,  -w1,    0,       0,    0,            0,
           0,       0,  Cab,    0,      w1, -k1b,    0,       0,    0,      M0b*R1b,
         Cac,       0,    0,    0,       0,    0, -k2c, -(wc-w),    0,            0,
           0,     Cac,    0,    0,       0,    0, wc-w,    -k2c,  -w1,            0,
           0,       0,  Cac,    0,       0,    0,    0,      w1, -k1c,      M0c*R1c;
    
    return (dM*duration).exp();
}



Eigen::Vector9d bm(
    Species const & species_a, Species const & species_b, Species const & species_c,
    double Cb, double Cc,
    double w0, double delta_w_rf, double w1, double duration,
    Eigen::Vector9d const & M)
{
    auto const M0a = species_a.M0;
    auto const R1a = 1/species_a.T1;
    auto const R2a = 1/species_a.T2;
    
    auto const M0b = species_b.M0;
    auto const R1b = 1/species_b.T1;
    auto const R2b = 1/species_b.T2;
    
    auto const M0c = species_c.M0;
    auto const R1c = 1/species_c.T1;
    auto const R2c = 1/species_c.T2;
    
    auto const Cab = M0b/M0a * Cb;
    auto const Cac = M0c/M0a * Cc;
    auto const Ca = Cab + Cac;
    
    auto const k1a = R1a+Ca;
    auto const k2a = R2a+Ca;
    
    auto const k1b = R1b+Cb;
    auto const k2b = R2b+Cb;
    
    auto const k1c = R1c+Cc;
    auto const k2c = R2c+Cc;
    
    auto const wa = w0*(1 + species_a.delta_w*1e-6);
    auto const wb = w0*(1 + species_b.delta_w*1e-6);
    auto const wc = w0*(1 + species_b.delta_w*1e-6);
    
    auto const w = w0*(1 + delta_w_rf*1e-6);
    
    Eigen::Matrix9d A;
    A << 
    //   Mxa      Mya   Mza   Mxb      Myb   Mzb   Mxc      Myc   Mzc
        -k2a, -(wa-w),    0,   Cb,       0,    0,   Cc,       0,    0,
        wa-w,    -k2a,  -w1,    0,      Cb,    0,    0,      Cc,    0,
           0,      w1, -k1a,    0,       0,   Cb,    0,       0,   Cc,
         Cab,       0,    0, -k2b, -(wb-w),    0,    0,       0,    0,
           0,     Cab,    0, wb-w,    -k2b,  -w1,    0,       0,    0,
           0,       0,  Cab,    0,      w1, -k1b,    0,       0,    0,
         Cac,       0,    0,    0,       0,    0, -k2c, -(wc-w),    0,
           0,     Cac,    0,    0,       0,    0, wc-w,    -k2c,  -w1,
           0,       0,  Cac,    0,       0,    0,    0,      w1, -k1c;
    
    Eigen::Vector9d const b{0, 0, M0a*R1a, 0, 0, M0b*R1b, 0, 0, M0c*R1c};
    
    Eigen::Vector9d const AinvB = A.partialPivLu().solve(b);
    
    return (A*duration).exp() * (M+AinvB) - AinvB;
}



Eigen::Vector10d
bm(
    Species const & species_a, Species const & species_b, Species const & species_c,
    double Cb, double Cc,
    double w0, double delta_w_rf, Eigen::VectorXd const & w1,
    double step, Eigen::Vector10d const & M0)
{
    Eigen::Vector10d M = M0;
    for(auto && w1_: w1)
    {
        M = bm(species_a, species_b, species_c, Cb, Cc, w0, delta_w_rf, w1_, step) * M;
    }
    return M;
}



Eigen::Vector9d
bm(
    Species const & species_a, Species const & species_b, Species const & species_c,
    double Cb, double Cc,
    double w0, double delta_w_rf, Eigen::VectorXd const & w1,
    double step, Eigen::Vector9d const & M0)
{
    Eigen::Vector9d M = M0;
    for(auto && w1_: w1)
    {
        M = bm(species_a, species_b, species_c, Cb, Cc, w0, delta_w_rf, w1_, step, M);
    }
    return M;
}


void bm3_partial(pybind11::module & m)
{
    using namespace pybind11::literals;
    
    m.def(
        "bm",
        pybind11::overload_cast<
            Species const &, Species const &, Species const &,
            double, double, double, double, double, double>(&bm),
        "species_a"_a, "species_b"_a, "species_c"_a,
        "Cb"_a, "Cc"_a, "w0"_a, "delta_w_rf"_a, "w1"_a, "step"_a);
    
    m.def(
        "bm",
        pybind11::overload_cast<
            Species const &, Species const &, Species const &,
            double, double, double, double, double, double, Eigen::Vector9d const &>(&bm),
        "species_a"_a, "species_b"_a, "species_c"_a, "Cb"_a, "Cc"_a, "w0"_a, "delta_w_rf"_a, "w1"_a,
        "step"_a, "M0"_a);
    
    m.def(
        "bm",
        pybind11::overload_cast<
            Species const &, Species const &, Species const &,
            double, double, double, double, Eigen::VectorXd const &, double,
            Eigen::Vector10d const &>(&bm),
        "species_a"_a, "species_b"_a, "species_c"_a, "Cb"_a, "Cc"_a, "w0"_a, "delta_w_rf"_a, "w1"_a,
        "step"_a, "M0"_a);
    
    m.def(
        "bm",
        pybind11::overload_cast<
            Species const &, Species const &, Species const &,
            double, double, double, double, Eigen::VectorXd const &, double,
            Eigen::Vector9d const &>(&bm),
        "species_a"_a, "species_b"_a, "species_c"_a, "Cb"_a, "Cc"_a, "w0"_a, "delta_w_rf"_a, "w1"_a,
        "step"_a, "M0"_a);
}
