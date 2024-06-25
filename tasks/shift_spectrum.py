import re

import nibabel
import numpy
import scipy
import spire

from . import utils

class ShiftSpectrum(spire.TaskFactory):
    def __init__(self, image, B0, shifted):
        spire.TaskFactory.__init__(self, str(shifted))
        self.file_dep = [image, B0]
        self.targets = [shifted]
        self.actions = [(__class__.action, (image, B0, shifted))]
    
    @staticmethod
    def action(image_path, B0_path, shifted_path):
        # Get the frequency information from the meta-data
        meta_data_path = re.sub("\.nii(\.gz)?$", ".json", str(image_path))
        ppm = utils.get_ppm(meta_data_path)
        
        # Load the image and the B0 map
        image, B0_image = [nibabel.load(x) for x in [image_path, B0_path]]
        data, B0_data = [numpy.array(x.dataobj) for x in [image, B0_image]]
        
        # Correct the nominal voxel-wise using the B0 map
        shifted_ppm = ppm[None, None, None, :]+B0_data[..., None]
        
        # Interpolate the data pixel-wise to shift the Z-spectrum
        shifted_ppm_flat = shifted_ppm.reshape(-1, shifted_ppm.shape[-1])
        data_flat = data.reshape(-1, data.shape[-1])
        shifted_data_flat = numpy.zeros_like(data_flat)
        for index, (x, fp) in enumerate(zip(shifted_ppm_flat, data_flat)):
            shifted_data_flat[index] = numpy.interp(x, ppm, fp)
        
        shifted_data = shifted_data_flat.reshape(data.shape)
        
        nibabel.save(
            nibabel.Nifti1Image(shifted_data, image.affine), shifted_path)
