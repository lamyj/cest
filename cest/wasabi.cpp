// WARNING: Stan must be included before Eigen so that the plugin system is
// active. https://discourse.mc-stan.org/t/includes-in-user-header/26093
#include <stan/math.hpp>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <xtensor-python/pyarray.hpp>

#include "slimp/actions.h"

#include "wasabi_sampler.h"

#include "univariate/predict_posterior.h"
#include "univariate/sampler.h"

using ArrayI = xt::xarray<int>;
using ArrayD = xt::xarray<double>;

std::tuple<ArrayD, std::vector<std::string>, std::vector<std::string>, ArrayI, ArrayI>
wasabi(
    ArrayD const & Delta_w, ArrayD const & Z,
    double B0, double B1_nominal, double t_p,
    slimp::action_parameters::Sample /* const & */ parameters)
{
    // Model type
    using Model = slimp::Model<wasabi_sampler::model>;
    
    // Model context, initialized with constant fields and dummy values for
    // non-constant fields.
    slimp::VarContext context;
    context.set("N", int(Z.shape()[1]));
    context.set("Delta_w", Delta_w);
    context.set("B0", B0);
    context.set("B1_nominal", B1_nominal);
    context.set("t_p", t_p);
    auto update_context = [&](slimp::VarContext & context, std::size_t r) {
        context.set("Z", xt::eval(xt::view(Z, r)));
    };
    update_context(context, 0);
    
    // Number of models to run
    std::size_t const R = Z.shape()[0];
    
    // Outputs: parameters and names, summarized HMC information
    ArrayD posterior;
    std::vector<std::string> names, raw_names;
    std::size_t depth_index, divergent_index;
    ArrayI exceeded_depth{ArrayI::shape_type{R}};
    ArrayI divergent{ArrayI::shape_type{R}};
    
    // Create a dummy model to initialize the result array and helpers
    {
        Model dummy(context, parameters);
        auto samples = dummy.create_samples();
        
        // Build full name vector to find the index of the treedepth and
        // divergent fields
        auto const hmc = dummy.hmc_names();
        depth_index =
            std::find(hmc.begin(), hmc.end(), "treedepth__") - hmc.begin();
        divergent_index = 
            std::find(hmc.begin(), hmc.end(), "divergent__") - hmc.begin();
        
        // We only store the model parameters, not the HMC parameters
        names = dummy.model_names();
        posterior = ArrayD(
            ArrayD::shape_type{R, names.size(), samples.shape(1), samples.shape(2)});
        
        raw_names = dummy.model_names(false, false);
    }
    
    auto update_results = [&](xt::xtensor<double, 3> const & samples, std::size_t r) {
        xt::view(posterior, r) = xt::view(
            samples, xt::range(-names.size(), xt::placeholders::_));
        
        exceeded_depth[r] = xt::sum(
                xt::view(samples, depth_index) > parameters.hmc.max_depth
            )[0];
        divergent[r] = xt::sum(xt::view(samples, divergent_index))[0];
    };
    
    stan::math::init_threadpool_tbb(parameters.threads_per_chain);
    
    slimp::parallel_sample<Model>(
        context, parameters, R, update_context, update_results);
    
    return {posterior, names, raw_names, exceeded_depth, divergent};
}

void wasabi(pybind11::module & m)
{
    m.def("wasabi", &slimp::sample<wasabi_sampler::model>);
    m.def(
        "wasabi",
        pybind11::overload_cast<
            ArrayD const &, ArrayD const &, double, double, double,
            slimp::action_parameters::Sample>(wasabi));
}
