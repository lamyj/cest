#include <Eigen/Core>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include "bm2.h"
#include "misc.h"

PYBIND11_MODULE(_cest, m)
{
    misc(m);
    bm2(m);
}
