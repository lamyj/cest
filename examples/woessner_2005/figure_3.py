import matplotlib.pyplot
import numpy

import cest

from common import Delta_omega_ppm, duration, M0, rf_offsets_Hz, rf_offsets_ppm, w0

Cb = 3000 # Hz
fb = 0.0007272
w1 = 2*numpy.pi * 128

T1s = [0.1, 0.3, 1, 2, 3, 5]
T2s = [0.003, 0.01, 0.03, 0.1, 0.3, 1]

species_b = cest.Species(0.1, 0.08, Delta_omega_ppm, fb)

bm = lambda rf_offset: cest.bm(
    species_a, species_b, Cb, w0, rf_offset, w1, duration, M)

data = numpy.empty((2, len(T1s), len(rf_offsets_ppm)))
for index, T1a in enumerate(T1s):
    species_a = cest.Species(T1a, 0.1, 0, 1)
    
    M = M0(species_a, species_b)
    data[0, index] = [bm(x)[2] for x in rf_offsets_ppm]
for index, T2a in enumerate(T2s):
    species_a = cest.Species(2, T2a, 0, 1)
    
    M = M0(species_a, species_b)
    data[1, index] = [bm(x)[2] for x in rf_offsets_ppm]

figure, plots = matplotlib.pyplot.subplots(
    2, 1, sharex=True, sharey=True, layout="constrained")
for d, plot, label, values in zip(data, plots, ["$T_1$", "$T_2$"], [T1s, T2s]):
    for i, (y, t) in enumerate(zip(d, values)):
        plot.plot(
            rf_offsets_Hz, y/species_a.M0, f"C{i}", lw=1,
            label=f"{label} = {t} s")
    plot.legend()

plots[1].set(
    xlim=(rf_offsets_Hz[0], rf_offsets_Hz[-1]), ylim=(0, 1),
    xlabel="RF frequency offset (Hz)")
plots[0].invert_xaxis()

figure.savefig("figure_3.png")
