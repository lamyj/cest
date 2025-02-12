import matplotlib.pyplot
import numpy

import cest

from common import Delta_omega_ppm, duration, M0, rf_offsets_Hz, rf_offsets_ppm, w0

fb = 0.0003636
w1 = 2*numpy.pi * 512

Cbs = [500, 5000]

species_a = cest.Species(1, 0.2, 0, 1)
species_b = cest.Species(0.1, 0.1, Delta_omega_ppm, fb)
M0 = M0(species_a, species_b)

bm = lambda Cb, rf_offset: cest.bm(
    species_a, species_b, Cb, w0, rf_offset, w1, duration, M0)

data = numpy.empty((len(Cbs), len(rf_offsets_ppm), 2))

for i_Cb, Cb in enumerate(Cbs):
    for i_rf_offset, rf_offset in enumerate(rf_offsets_ppm):
        M = bm(Cb, rf_offset)
        data[i_Cb, i_rf_offset, :] = M[2], M[5]

figure, plots = matplotlib.pyplot.subplots(
    2, 1, sharex=True, sharey=True, layout="constrained")
for i, (plot, d, Cb) in enumerate(zip(plots, data, Cbs)):
    plot.plot(rf_offsets_Hz, d[:, 0]/species_a.M0, "black", lw=1.5)
    plot.plot(rf_offsets_Hz, d[:, 1]/species_b.M0, "black", lw=1)
    plot.text(15000, 0.5, f"$C_b$ = {Cb} s$^{{-1}}$")

plots[0].set(xlim=(rf_offsets_Hz[0], rf_offsets_Hz[-1]), ylim=(0, 1))
plots[0].invert_xaxis()
plots[1].set(xlabel="RF frequency offset (Hz)")

figure.savefig("figure_1.png")
