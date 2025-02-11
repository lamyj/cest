from ._cest import Species

water = Species(1.2, 40e-3, 0, 1, 0)

# Values from "Analysis of chemical exchange saturation transfer contributions
# from brain metabolites to the Z-spectra at various field strengths and pH",
# Khlebnikov et al. Scientific Reports 9(1), 2019.
# doi:10.1038/s41598-018-37295-y
# Notes
# - concentrations in the paper mentioned above are given in mM. The 1/55 factor
#   comes from the molarity of pure water
# - the paper mentions that the T2 of labile protons of taurine and mobile
#   amides could not be fitted reliably and was set at 10 ms
myo_inositol = Species(1.2, 22.8e-3, 1, 0.0054/55, 2090)
creatine = Species(1.2, 7.1e-3, 2, 0.00705/55, 810)
phosphocreatine_1_93ppm = Species(1.2, 7.8e-3, 1.93, 0.00705/55, 67)
phosphocreatine_2_64ppm = Species(1.2, 7.8e-3, 2.64, 0.00705/55, 126)
gaba = Species(1.2, 17.2e-3, 2.91, 0.0015/55, 6900)
taurine = Species(1.2, 10e-3, 3.18, 0.00155/55, 49600)
glutamate = Species(1.2, 6.9e-3, 3.2, 0.0066/55, 7480)
glutamine_2_15ppm = Species(1.2, 13.8e-3, 2.15, 0.003/55, 17)
glutamine_2_87ppm = Species(1.2, 13.8e-3, 2.87, 0.003/55, 49)
glutamine_3_18ppm = Species(1.2, 13.8e-3, 3.18, 0.003/55, 22880)
mobile_amides = Species(1.2, 10e-3, 3.5, 0.072/55, 22)
