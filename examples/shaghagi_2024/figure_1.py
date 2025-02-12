import matplotlib.pyplot
import numpy

import cest

gamma = 267522187.44
gamma_bar = gamma/(2*numpy.pi)

B0 = 3 # T
w0 = gamma * B0 # rad/s

# Water is 55 mol/L, and has two protons
species_a = cest.Species(2, 0.1, 0, 2*55)
species_b = cest.Species(1, 0.1, 3.5, 0.1)

Cb = 30 # Hz

B1s = 1e-6*numpy.array([0.5, 1, 1.5]) # T
w1s = gamma * B1s # rad/s

delta_w_rfs = numpy.arange(-10, 10+0.125, 0.125)
duration = 10 # s

spectra = numpy.empty((len(w1s), len(delta_w_rfs)))
for index, w1 in enumerate(w1s):
    spectra[index] = [
        cest.bm(
                species_a, species_b, Cb, w0, delta_w_rf, w1, duration,
                [0, 0, species_a.M0, 0, 0, species_b.M0]
            )[2]
        for delta_w_rf in delta_w_rfs]

figure, plot = matplotlib.pyplot.subplots(layout="constrained")
for spectrum, B1, color in zip(spectra, B1s, ["green", "red", "blue"]):
    plot.plot(
        delta_w_rfs, spectrum/species_a.M0, color, lw=1,
        label=f"$B_1$={B1*1e6} µT")

plot.set(
    xlim=(delta_w_rfs[0], delta_w_rfs[-1]), ylim=(0, 1),
    xlabel="$\Delta\omega_{RF}$ (ppm)", ylabel="Z")
plot.invert_xaxis()
plot.legend()

figure.savefig("figure_1.png")
