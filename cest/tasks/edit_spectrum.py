import json

import nibabel
import numpy
import spire

from . import utils

@spire.task_factory
def EditSpectrum(
        source_image: spire.file_dep, source_meta_data: spire.file_dep,
        operations,
        target_image: spire.target, target_meta_data: spire.target):
    source_image = nibabel.load(source_image)
    data = numpy.asarray(source_image.dataobj)
    ppm = utils.get_ppm(source_meta_data)
    
    for name, *operands in operations:
        if name == "range":
            (low, high), inclusion = operands
            inclusion = inclusion.lower() == "include"
            
            selector = (ppm >= low) & (ppm <= high)
            if not inclusion:
                selector = ~selector
            
            ppm = ppm[selector]
            data = data[..., selector]
        else:
            raise Exception(f"Unknown operation {name!r}")
    
    nibabel.save(nibabel.Nifti1Image(data, source_image.affine), target_image)
    
    with open(target_meta_data, "w") as fd:
        json.dump({"SaturationPulse": [{"FrequencyOffset": x} for x in ppm]}, fd)
