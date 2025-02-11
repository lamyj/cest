import matplotlib.pyplot
import numpy

import cest

# Experimental conditions: gyromagnetic ratio (rad/s/T), main magnetic field (Hz)
gamma = 267522187.44
B0 = 300e6

# CW saturation pulse: amplitude (T), duration (s)
B1_cw = 6e-6
tau_cw = 1

# Simulation parameters: time step (s), offset of the RF pulse (PPM), main species
step = 1e-3 # s
w_ppm = numpy.linspace(-5, 5, 101)
species_a = cest.species.water

# Shaped pulse. The RF pulses are normalized: multiply by γ⋅B₁⋅steps to get a
# saturation equivalent to the CW pulse defined above.
steps = int(round(tau_cw/step)) # Unitless
pulse = cest.pulses.gaussian(steps) * gamma*B1_cw * steps

# Compute the Z spectra of three species exchanging with water when applying the
# shaped pulse at the specified offsets
z_spectra = {}
for name in ["glutamate", "creatine", "mobile_amides"]:
    species_b = getattr(cest.species, name)
    # WARNING: unrealistic concentration, to better show the effect
    species_b.M0 *= 20
    
    z_spectra[name] = numpy.zeros_like(w_ppm)
    
    # NOTE: the frequency offset must be converted to Hz
    for i, w in enumerate(w_ppm*1e-6 * B0):
        M0 = numpy.array([0, 0, species_a.M0, 0, 0, species_b.M0, 1])
        M = cest.two_pools(species_a, species_b, w, pulse, step, B0, M0)
        z_spectra[name][i] = M[2]

mtrs = {}
for name, z_spectrum in z_spectra.items():
    mtrs[name] = cest.mtr(z_spectrum, w_ppm)
    
labels = {"glutamate": "Glu", "creatine": "Cr", "mobile_amides": "mAmides"}

figure, plot = matplotlib.pyplot.subplots(layout="constrained")
for name, z_spectrum in z_spectra.items():
    plot.plot(w_ppm, z_spectrum, label=labels[name])

plot.set(ylim=0, xlabel=r"$\Delta\omega$ (ppm)", ylabel="$M_{z,a}$ (unitless)")
plot.legend()
plot.invert_xaxis()

figure.savefig("z_spectra.png")

figure, plot = matplotlib.pyplot.subplots(layout="constrained")
for name, mtr in mtrs.items():
    plot.plot(w_ppm[w_ppm>=0], mtr, label=labels[name])

plot.set(ylim=0, xlabel=r"$\Delta\omega$ (ppm)", ylabel="MTR (unitless)")
plot.legend()
plot.invert_xaxis()

figure.savefig("mtr.png")
