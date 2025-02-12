import matplotlib.pyplot
import numpy

import cest

gamma = 267522187.44
gamma_bar = gamma/(2*numpy.pi)

B0 = 9.4 # T
w0 = gamma * B0 # rad/s

# Water is 55 mol/L, and has two protons
species_a = cest.Species(2, 0.1, 0, 2*55)
species_b = cest.Species(1, 0.1, 0, 1)

B1 = 1e-6 # T
w1 = gamma * B1 # rad/s
delta_w_rfs = numpy.linspace(-5, +5, 201)
duration = 10

#########################################
# Figure 2A: effect of Δω for species B #
#########################################

delta_w_bs = numpy.array([1, 2, 3.5])
species_b.M0 = 0.1
Cb = 30

delta_w_spectra = numpy.empty((len(delta_w_bs), len(delta_w_rfs)))
for i, delta_w_b in enumerate(delta_w_bs):
    species_b.delta_w = delta_w_b
    delta_w_spectra[i,:] = [
        cest.bm(
               species_a, species_b, Cb, w0, delta_w_rf, w1, duration,
               [0, 0, species_a.M0, 0, 0, species_b.M0]
            )[2]
        for delta_w_rf in delta_w_rfs]

###########################
# Figure 2B: effect of Cb #
###########################

species_b.delta_w=3.5
species_b.M0 = 0.1
Cbs = numpy.array([50, 500, 5000])

Cb_spectra = numpy.empty((len(Cbs), len(delta_w_rfs)))
for i, Cb in enumerate(Cbs):
    Cb_spectra[i,:] = [
        cest.bm(
               species_a, species_b, Cb, w0, delta_w_rf, w1, duration,
               [0, 0, species_a.M0, 0, 0, species_b.M0]
            )[2]
        for delta_w_rf in delta_w_rfs]

#########################################
# Figure 2C: effect of M0 for species B #
#########################################

species_a.M0 = 110
species_b.delta_w=3.5
M0s = numpy.array([1e-3, 10e-3, 100e-3])
Cb = 500

M0_spectra = numpy.empty((len(M0s), len(delta_w_rfs)))
for i, M0 in enumerate(M0s):
    species_b.M0 = M0
    M0_spectra[i,:] = [
        cest.bm(
               species_a, species_b, Cb, w0, delta_w_rf, w1, duration,
               [0, 0, species_a.M0, 0, 0, species_b.M0]
            )[2]
        for delta_w_rf in delta_w_rfs]

##########################
# Plot the three effects #
##########################

figure, plots = matplotlib.pyplot.subplots(
    1, 3, sharex=True, sharey=True, layout="constrained", figsize=(9, 3))

colors = ["green", "red", "blue"]
for spectrum, delta_w, color in zip(delta_w_spectra, delta_w_bs, colors):
    plots[0].plot(
        delta_w_rfs, 1-spectrum/species_a.M0, color,
        label=f"$\Delta\omega$={delta_w} ppm")
for spectrum, Cb, color in zip(Cb_spectra, Cbs, colors):
    plots[1].plot(
        delta_w_rfs, 1-spectrum/species_a.M0, color,
        label=f"Cb={Cb} $s^{{-1}}$")
for spectrum, M0, color in zip(M0_spectra, M0s, colors):
    plots[2].plot(
        delta_w_rfs, 1-spectrum/species_a.M0, color, label=f"$M_0$={M0} mmol")

plots[0].set(
    xlim=(delta_w_rfs[0], delta_w_rfs[-1]), ylim=(0, 1), ylabel="1-Z")
for plot in plots:
    plot.set(xticks=range(-5, 6), xlabel="$\Delta\omega_{RF}$ (ppm)")
    plot.legend(prop={"size": 7})

figure.savefig("figure_2.png")
