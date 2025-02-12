import matplotlib.pyplot
import numpy

import cest

from common import Delta_omega_ppm, duration, M0, rf_offsets_Hz, rf_offsets_ppm, w0

fb = 0.0003636
w1 = 2*numpy.pi * 512

Cbs = [300, 500, 900, 3000, 7000, 9000]

species_a = cest.Species(2, 0.2, 0, 1)
species_b = cest.Species(0.1, 0.1, Delta_omega_ppm, fb)
M0 = M0(species_a, species_b)

bm = lambda Cb, rf_offset: cest.bm(
    species_a, species_b, Cb, w0, rf_offset, w1, duration, M0)

data = numpy.empty((len(Cbs), len(rf_offsets_ppm)))

for i_Cb, Cb in enumerate(Cbs):
    for i_rf_offset, rf_offset in enumerate(rf_offsets_ppm):
        data[i_Cb, i_rf_offset] = bm(Cb, rf_offset)[2]

figure, plot = matplotlib.pyplot.subplots(
    sharex=True, sharey=True, layout="constrained")
for i, (d, Cb) in enumerate(zip(data, Cbs)):
    plot.plot(
        rf_offsets_Hz, d/species_a.M0, f"C{i}", lw=1,
        label=f"$C_b$ = {Cb} s$^{{-1}}$")

plot.set(
    xlim=(rf_offsets_Hz[0], rf_offsets_Hz[-1]), ylim=(0, 1),
    xlabel="RF frequency offset (Hz)")
plot.invert_xaxis()
plot.legend()

figure.savefig("figure_2.png")
