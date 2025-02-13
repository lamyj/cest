Simulation
==========

Simulation of CEST is based on the Bloch-McConnell equation. You will first need to define which species you want to simulate. Each species is a instance of the :py:class:`Species <cest.Species>` class.

A single saturation block-pulse can be simulated as follows, with a set of frequency offsets to obtain a Z-spectrum. This examples is based on the two-pools Bloch-McConnell equations, using the formalism of `Numerical solution of the Bloch equations provides insights into the optimum design of PARACEST agents for MRI <https://doi.org/10.1002/mrm.20408>`_ (Woessner et al. Magnetic Resonance in Medicine, 53(4), 2005).

.. code-block:: python
    
    import cest
    import matplotlib.pyplot
    import numpy
    
    # Experimental conditions: gyromagnetic ratio (rad/s/T), main magnetic field (T and rad/s)
    gamma = 267522187.44 # rad/s/T
    B0 = 3 # T
    w0 = gamma * B0 # rad/s
    
    # Block saturation pulse: amplitude and duration
    B1 = 12e-6 # T
    w1 = gamma * B1 # rad/m
    tau = 10 # s
    
    # Species and their exchange rate (B -> A)
    species_a = cest.Species(1, 0.2, 0, 1)
    species_b = cest.Species(0.1, 0.1, 160, 0.0003636)
    Cb = 500 # Hz
    
    # Simulate between -100 and +200 ppms
    offsets = numpy.linspace(-100, +200, 401)
    magnetization = numpy.array([
        cest.bm(
            species_a, species_b, Cb, w0, offset, w1, tau,
            [0, 0, species_a.M0, 0, 0, species_b.M0])
        for offset in offsets])
    
    # Plot the results
    matplotlib.pyplot.plot(offsets, magnetization[:,2]/species_a.M0, lw=1)
    matplotlib.pyplot.xlabel("$\Delta\omega_{RF}$ (ppm)")
    matplotlib.pyplot.ylim(0, 1)
    matplotlib.pyplot.ylabel("$M_z/M_0$")
    matplotlib.pyplot.gca().invert_xaxis()

Shaped pulses can be defined using the `pre-defined shapes <api/functions.html#pulses>`__ or by adding your own. Each shape is normalized and discretized, so it needs to be scaled:

.. code-block:: python
    
    # Discretization step
    step = 10e-3 # s
    
    # Shaped pulse. The RF pulses are normalized: multiply by ω₁⋅N to get a
    # saturation equivalent to the block pulse defined above.
    steps = int(round(tau/step)) # Unitless
    pulse = cest.pulses.gaussian(steps) * w1 * steps # rad/s

The shaped pulse is then simulated using the same function as above:

.. code-block:: python
    
    magnetization = numpy.array([
        cest.bm(
            species_a, species_b, Cb, w0, offset, pulse, step,
            [0, 0, species_a.M0, 0, 0, species_b.M0])
        for offset in offsets])
    
    # Plot the results
    matplotlib.pyplot.plot(offsets, magnetization[:,2]/species_a.M0, lw=1)
