#include <cmath>
#include <iostream>

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace Eigen
{
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
    double w;
    /// @brief Equilibrium magnetization (unitless)
    double M0;
    /// @brief Exchange rate in Hz
    double k;
};

/**
 * @brief Simulation using the two pools Bloch-McConnel model
 * @param w Frequency offset of the saturation RF pulse in Hz
 * @param w1 Frequency of the B1 field of the saturation RF pulse in Hz
 * @param duration Duration of the simulation in s
 * @param B0 Nominal B0 field in Hz
 * @param Delta_B0 B0 offset in PPM
 */
Eigen::Matrix7d
two_pools(
    Species const & species_a, Species const & species_b,
    double w, double w1, double duration, double B0, double Delta_B0=0)
{
    auto const R1a = 1/species_a.T1;
    auto const R2a = 1/species_a.T2;
    auto const R1b = 1/species_b.T1;
    auto const R2b = 1/species_b.T2;
    
    auto const wa = species_a.w*1e-6 * B0;
    auto const wb = species_b.w*1e-6 * B0;
    
    auto const M0a = species_a.M0;
    auto const M0b = species_b.M0;
    
    auto const fb = species_b.M0/species_a.M0;
    
    auto const ka = species_b.k;
    
    auto const w0 = Delta_B0*1e-6 * B0;
    
    auto const omega_a = 2*M_PI * (w-wa+w0);
    auto const omega_b = 2*M_PI * (w-wb+w0);
    
    auto const fbka = fb*ka;
    
    Eigen::Matrix7d A;
    A << 
        -R2a-fbka,  -omega_a,         0,      ka,        0,       0,       0,
          omega_a, -R2a-fbka,       -w1,       0,       ka,       0,       0,
                0,        w1, -R1a-fbka,       0,        0,      ka, R1a*M0a,
             fbka,         0,         0, -R2b-ka, -omega_b,       0,       0,
                0,      fbka,         0, omega_b,  -R2b-ka,     -w1,       0,
                0,         0,      fbka,       0,       w1, -R1b-ka, R1b*M0b,
                0,         0,         0,       0,        0,       0,       0;
    
    return (A*duration).exp();
}

Eigen::Vector7d
two_pools(
    Species const & species_a, Species const & species_b,
    double w, Eigen::VectorXd const & w1, double step, double B0,
    Eigen::Vector7d const & M0, double Delta_B0=0)
{
    Eigen::Vector7d M = M0;
    for(auto && w1_: w1)
    {
        M = two_pools(species_a, species_b, w, w1_, step, B0, Delta_B0) * M;
    }
    return M;
}

PYBIND11_MODULE(_cest, m)
{
    using namespace pybind11::literals;
    
    pybind11::class_<Species>(m, "Species", "Chemical species")
        .def_readwrite("T1", &Species::T1, "T1 relaxation times in s")
        .def_readwrite("T2", &Species::T2, "T2 relaxation times in s")
        .def_readwrite("w", &Species::w, "Frequency offset in ppm")
        .def_readwrite("M0", &Species::M0, "Equilibrium magnetization (unitless)")
        .def_readwrite("k", &Species::k, "Exchange rate in Hz")
        .def(
            pybind11::init(
                [](double T1, double T2, double w, double M0, double k){
                    return Species{T1, T2, w, M0, k};
                }),
            "T1"_a, "T2"_a, "w"_a, "M0"_a, "k"_a);
    m.def(
        "two_pools",
        pybind11::overload_cast<
                Species const &, Species const &, double, 
                double,
                double, double,
                double>(&two_pools),
        "species_a"_a, "species_b"_a, "w"_a,
        "w1"_a,
        "duration"_a, "B0"_a,
        "delta_B0"_a=0);
    m.def(
        "two_pools",
        pybind11::overload_cast<
                Species const &, Species const &, double,
                Eigen::VectorXd const &,
                double, double,
                Eigen::Vector7d const &, double>(&two_pools),
        "species_a"_a, "species_b"_a, "w"_a,
        "w1"_a,
        "duration"_a, "B0"_a,
        "M0"_a, "delta_B0"_a=0);
}
