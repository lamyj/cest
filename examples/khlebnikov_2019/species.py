import dataclasses

import cest

# Water is 55 mol/L, and has two protons: M0 is 110

white_matter = {
    3: cest.Species(0.9, 70e-3, 0, 110),
    4.7: cest.Species(1.07, 60e-3, 0, 110),
    7: cest.Species(1.2, 40e-3, 0, 110),
    9.4: cest.Species(1.66, 37, 0, 110),
    11.7: cest.Species(1.74, 27e-3, 0, 110),
    14.1: cest.Species(1.98, 20.4e-3, 0, 110)
}

grey_matter = {
    3: cest.Species(1.45, 90e-3, 0, 110),
    4.7: cest.Species(1.47, 73e-3, 0, 110),
    7: cest.Species(1.8, 55e-3, 0, 110),
    9.4: cest.Species(2.1, 42e-3, 0, 110),
    11.7: cest.Species(2.11, 37e-3, 0, 110),
    14.1: cest.Species(2.45, 25.9e-3, 0, 110)
}

@dataclasses.dataclass
class LabileProton:
    delta_w: float
    Cb: float

@dataclasses.dataclass
class Metabolite:
    T2: float
    labile_protons: list[LabileProton]
    concentration: float
    
    def species(self, T1, proton_index):
        p = self.labile_protons[proton_index]
        return cest.Species(T1, self.T2, p.delta_w, self.concentration)

Glc = Metabolite(
    6.9e-3, [LabileProton(*x) for x in [
        (2.18, 3860), (2.88, 4750), (1.10, 3940), (1.39, 2560), (0.74, 950)]],
    1e-3)
MI = Metabolite(22.8e-3, [LabileProton(1.0, 2090)], 5.9e-3)
Cr = Metabolite(7.1e-3, [LabileProton(2, 810)], 8.4e-3)
PCr = Metabolite(7.8e-3, [LabileProton(2.64, 126)], 8.4e-3)
GABA = Metabolite(17.2e-3, [LabileProton(2.91, 6900)], 1.5e-3)
Tau = Metabolite(10e-3, [LabileProton(3.18, 49600)], 2.1e-3)
Glu = Metabolite(6.9e-3, [LabileProton(3.2, 7480)], 11.5e-3)
Gln = Metabolite(
    13.8e-3, [LabileProton(*x) for x in [
        (2.15, 17), (2.87, 49), (3.18, 22880)]],
    3e-3)
NAA = Metabolite(None, [LabileProton(3.33, 0)], 12e-3)
mAmides = Metabolite(10e-3, [LabileProton(3.5, 22)], 72e-3)
