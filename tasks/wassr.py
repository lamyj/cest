import re

import nibabel
import numpy
import scipy
import spire

from . import utils

class WASSR(spire.TaskFactory):
    def __init__(self, image, wassr, delta_ppm=0.001, ppm_range=None):
        spire.TaskFactory.__init__(self, str(wassr))
        self.file_dep = [image]
        self.targets = [wassr]
        self.actions = [(__class__.action, (image, wassr, delta_ppm, ppm_range))]
    
    @staticmethod
    def action(image_path, wassr_path, delta_ppm=0.001, ppm_range=None):
        # Get the frequency information from the meta-data
        meta_data_path = re.sub("\.nii(\.gz)?$", ".json", str(image_path))
        ppm = utils.get_ppm(meta_data_path)
        
        # Load the image data
        image = nibabel.load(image_path)
        data = numpy.array(image.dataobj)
        
        # Sort the volumes in increasing PPM order
        order = numpy.argsort(ppm)
        ppm = ppm[order]
        data = data[..., order]
        
        # Keep the requested data
        if ppm_range:
            selector = (ppm >= ppm_range[0]) & (ppm <= ppm_range[1])
            ppm = ppm[selector]
            data = data[..., selector]
        
        # Interpolate the Z-spectrum
        ppm_fine = numpy.arange(ppm.min(), ppm.max()+delta_ppm, delta_ppm)
        interpolator = scipy.interpolate.CubicSpline(ppm, data, axis=3)
        interpolated = interpolator(ppm_fine)
        
        # Compute the B0 map from the minimum of the interpolated Z-spectrum
        min_index = numpy.argmin(interpolated, axis=3)
        B0_ppm = ppm_fine[min_index]

        nibabel.save(nibabel.Nifti1Image(B0_ppm, image.affine), wassr_path)



