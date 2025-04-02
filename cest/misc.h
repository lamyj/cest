#ifndef _94ad7880_73fb_49c3_87f2_2f780574a299
#define _94ad7880_73fb_49c3_87f2_2f780574a299

#include <Eigen/Core>

#include <pybind11/pybind11.h>

namespace Eigen
{
using Matrix6d = Eigen::Matrix<double, 6, 6>;
using Vector6d = Eigen::Matrix<double, 6, 1>;
using RowVector6d = Eigen::Matrix<double, 1, 6>;

using Matrix7d = Eigen::Matrix<double, 7, 7>;
using Vector7d = Eigen::Matrix<double, 7, 1>;
using RowVector7d = Eigen::Matrix<double, 1, 7>;

using Matrix9d = Eigen::Matrix<double, 9, 9>;
using Vector9d = Eigen::Matrix<double, 9, 1>;
using RowVector9d = Eigen::Matrix<double, 1, 9>;

using Matrix10d = Eigen::Matrix<double, 10, 10>;
using Vector10d = Eigen::Matrix<double, 10, 1>;
using RowVector10d = Eigen::Matrix<double, 1, 10>;
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

void misc(pybind11::module m);

#endif // _94ad7880_73fb_49c3_87f2_2f780574a299
