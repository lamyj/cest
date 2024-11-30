Simulation
==========

Simulation of CEST is based on the Bloch-McConnell equation. You will first need to define which species you want to simulate. Each species is a instance of the :py:class:`Species <cest.Species>` class, with a few pre-defined ones (e.g. `cest.species.water`, `cest.species.glutamate`).

You will then need to define a saturation pulse. Shapes can be defined using the `pre-defined shapes <api/functions.html#pulses>`__ or by adding your own. Each shape is normalized and discretized, so it needs to be scaled:

.. code-block:: python
    
    import cest
    import numpy
    
    # Experimental conditions: gyromagnetic ratio (rad/s/T), main magnetic field (Hz)
    gamma = 267522187.44
    B0 = 300e6
    
    # CW saturation pulse: amplitude (T), duration (s)
    B1_cw = 6e-6
    tau_cw = 1
    
    # Simulation time step (s)
    step = 1e-3
    
    # Shaped pulse. The RF pulses are normalized: multiply by γ⋅B₁⋅N to get a
    # saturation equivalent to the CW pulse defined above.
    steps = int(round(tau_cw/step)) # Unitless
    pulse = cest.pulses.gaussian(steps) * gamma*B1_cw * steps

Individual pulse steps can the be simulated with the :py:func:`two_pools <cest.two_pools>` function:

.. code-block:: python
    
    # Offset of the saturation RF pulse (PPM)
    w_ppm = -3.5
    
    # Magnetization of the two pools
    M = numpy.array([
        0, 0, cest.species.water.M0,
        0, 0, cest.species.glutamate.M0,
        1])
    for w1 in pulse:
        M = cest.two_pools(
                cest.species.water, cest.species.glutamate,
                w_ppm*1e-6 * B0, w1, step, B0
            ) @ M
    print(M[2])
