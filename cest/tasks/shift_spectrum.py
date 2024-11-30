import re

import nibabel
import numpy
import scipy
import spire

from . import utils

class ShiftSpectrum(spire.TaskFactory):
    """Shift a Z-spectrum as described in the WASSR method.
    
    Parameters
    ----------
    
    image : path_like
        Path to the source Z-spectrum image
    meta_data : path_like
        Path to the meta-data related to the source image
    B0 : path_like
        Path to the frequency shift map
    shifted : path_like
        Path to the target shifted Z-spectrum
    
    References
    ----------
    
    *Water saturation shift referencing (WASSR) for chemical \
        exchange saturation transfer (CEST) experiments*, Kim et al., Magnetic \
        Resonance in Medicine 61(6), 2009. \
        `doi:10.1002/mrm.21873 <https://doi.org/10.1002/mrm.21873>`_.
    """
    
    def __init__(self, image, meta_data, B0, shifted):
        spire.TaskFactory.__init__(self, str(shifted))
        self.file_dep = [image, meta_data, B0]
        self.targets = [shifted]
        self.actions = [(__class__.action, (image, meta_data, B0, shifted))]
    
    @staticmethod
    def action(image, meta_data, B0, shifted):
        # Get the frequency information from the meta-data
        ppm = utils.get_ppm(meta_data)
        
        # Load the image and the B0 map
        image, B0_image = [nibabel.load(x) for x in [image, B0]]
        data, B0_data = [numpy.array(x.dataobj) for x in [image, B0_image]]
        
        # Correct the nominal voxel-wise using the B0 map
        shifted_ppm = ppm[None, None, None, :]+B0_data[..., None]
        
        # Interpolate the data pixel-wise to shift the Z-spectrum
        shifted_ppm_flat = shifted_ppm.reshape(-1, shifted_ppm.shape[-1])
        data_flat = data.reshape(-1, data.shape[-1])
        shifted_data_flat = numpy.zeros_like(data_flat)
        # NOTE: is linear interpolation the best choice?
        for index, (x, fp) in enumerate(zip(shifted_ppm_flat, data_flat)):
            shifted_data_flat[index] = numpy.interp(x, ppm, fp)
        
        shifted_data = shifted_data_flat.reshape(data.shape)
        
        nibabel.save(nibabel.Nifti1Image(shifted_data, image.affine), shifted)
