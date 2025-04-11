#include <Eigen/Core>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#define FORCE_IMPORT_ARRAY
#include <xtensor-python/pyarray.hpp>

#include "bm2.h"
#include "bm3_partial.h"
#include "misc.h"
#include "wasabi.h"

PYBIND11_MODULE(_cest, m)
{
    misc(m);
    bm2(m);
    bm3_partial(m);
    wasabi(m);
}
