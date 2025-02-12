import itertools

import colorcet
import matplotlib.pyplot
import numpy

import cest

import species

gamma = 267522187.44
gamma_bar = gamma/(2*numpy.pi)

B0 = 7 # T
w0 = gamma * B0 # rad/s

species_a = species.grey_matter[B0]

metabolite_names = ["MI", "Cr", "PCr", "Glu", "mAmides"]
metabolites = [(getattr(species, x), 0) for x in metabolite_names]

B1s = numpy.linspace(0.1e-6, 10e-6, 100) # T
w1s = gamma * B1s # rad/s
durations = numpy.linspace(0.1, 5, 100) # s
mtrs = numpy.zeros((len(metabolites), len(w1s), len(durations)))

species_a = species.grey_matter[B0]

iterator = zip(
    numpy.ndindex(mtrs.shape), itertools.product(metabolites, w1s, durations))
for index, ((metabolite, proton), w1, duration) in iterator:
    species_b = metabolite.species(species_a.T1, proton)
    
    bm = lambda delta_w, w1, duration: cest.bm(
        species_a, species_b, metabolite.labile_protons[proton].Cb, w0,
        delta_w, w1, duration, [0, 0, species_a.M0, 0, 0, species_b.M0])[2]
    
    mtrs[index] += (
        (
            bm(-species_b.delta_w, w1, duration)
            - bm(+species_b.delta_w, w1, duration))
        / bm(100, w1, duration))

figure, plots = matplotlib.pyplot.subplots(
    3, 2, layout="constrained", figsize=(9, 8))
x, y = numpy.meshgrid(B1s*1e6, durations)
colormap = matplotlib.cm.get_cmap("viridis")
used_plots = [plots[0,0], plots[1,0], plots[1,1], plots[2,0], plots[2,1]]

xticks=[0.1, 2, 4, 6, 8, 10]
yticks=[0.1, 1, 2, 3, 4, 5]
for mtrs_, metabolite_name, plot in zip(mtrs, metabolite_names, used_plots):
    plot.contourf(x, y, mtrs_.T, 10, cmap=colormap)
    plot.contour(x, y, mtrs_.T, 10, colors="black", linewidths=1)
    
    optimum_index = numpy.unravel_index(numpy.argmax(mtrs_), mtrs_.shape)
    optimum = 1e6*B1s[optimum_index[0]], durations[optimum_index[1]]
    
    plot.plot(*optimum, "k+")
    plot.set(
        xticks=xticks, yticks=yticks,
        title=f"{metabolite_name} ({optimum[0]:.1f} µT, {optimum[1]:.1f} s)")
    labels = plot.get_xticklabels(), plot.get_yticklabels()
    plot.set(xticklabels=[], yticklabels=[])

for plot in plots[-1]:
    plot.set(xticklabels=xticks, xlabel="$B_1$ (µt)")
for plot in plots[:,0]:
    plot.set(yticklabels=yticks, ylabel="Saturation time (s)")

figure.colorbar(
    matplotlib.cm.ScalarMappable(
        norm=matplotlib.colors.Normalize(0, 1), cmap=colormap),
    cax=plots[0,1].inset_axes([0, 0.45, 1, 0.1]), orientation="horizontal",
    label="$MTR_{asym}$ (unitless)")
plots[0,1].set(xticks=[], yticks=[])
plots[0,1].spines[:].set_visible(False)

figure.savefig("figure_1.png")
