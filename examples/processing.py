import re

import numpy
import pathlib
import spire.ants
import yaml

import cest.tasks

root = pathlib.Path(".")

with open(root/"series.yml") as fd:
    series = yaml.load(fd, Loader=yaml.SafeLoader)

original = root/"original"

images = {n: original/s/"1.nii.gz" for n, s in series.items()}
meta_data = {n: re.sub("\.nii(\.gz)?", ".json", str(i)) for n, i in images.items()}

if "Water" in series:
    derived = root/"derived"/"WASSR"
    derived.mkdir(parents=True, exist_ok=True)

    B0 = cest.tasks.WASSR(images["Water"], meta_data["Water"], derived/"B0.nii.gz")
    B0_in_Glutamate = spire.ants.ApplyTransforms(
        B0.targets[0], (images["Glutamate"], 0), [],
        derived/"B0_in_Glutamate.nii.gz")
    shifted = cest.tasks.ShiftSpectrum(
        images["Glutamate"], meta_data["Glutamate"], B0_in_Glutamate.targets[0],
        derived/"Glutamate_shifted.nii.gz")
    
    ppms = numpy.linspace(-5, 5, 501)
    refined = cest.tasks.Refine(
        shifted.targets[0], shifted.file_dep[1], numpy.linspace(-5, 5, 501),
        derived/"Glutamate_refined.nii.gz")

derived = root/"derived"/"dummy"
derived.mkdir(parents=True, exist_ok=True)

B0 = cest.tasks.WASSR(
    images["Glutamate"], meta_data["Glutamate"], derived/"B0.nii.gz",
    ppm_range=(-5.1, +5.1))
shifted = cest.tasks.ShiftSpectrum(
    images["Glutamate"], meta_data["Glutamate"], B0.targets[0],
    derived/"Glutamate_shifted.nii.gz")

ppms = numpy.linspace(-5, 5, 501)
refined = cest.tasks.Refine(
    shifted.targets[0], shifted.file_dep[1], ppms,
    derived/"Glutamate_refined.nii.gz")

mtr = cest.tasks.MTR(refined.targets[0], ppms, derived/"Glutatmate_MTR.nii.gz")
