import re

import nibabel
import numpy
import scipy
import spire

from . import utils

class WASSR(spire.TaskFactory):
    """Compute a Z-spectrum shift map using the [WASSR]_ method.
    
    Parameters
    ----------
    image : path_like
        Path to the source image
    meta_data : path_like
        Path to the meta-data related to the source image
    wassr : path_like
        Path to the target WASSR map
    delta_ppm : float, optional
        Interpolation step
    ppm_range : pair_of_floats, optional
        Range over which to interpolate the Z-spectrum, defaults to full range
        defined by meta-data.
    
    References
    ----------
    .. [WASSR] *Water saturation shift referencing (WASSR) for chemical \
        exchange saturation transfer (CEST) experiments*, Kim et al., Magnetic \
        Resonance in Medicine 61(6), 2009. \
        `doi:10.1002/mrm.21873 <https://doi.org/10.1002/mrm.21873>`_.
    
    """
        
    def __init__(self, image, meta_data, wassr, delta_ppm=0.001, ppm_range=None):
        spire.TaskFactory.__init__(self, str(wassr))
        self.file_dep = [image, meta_data]
        self.targets = [wassr]
        self.actions = [
            (__class__.action, (image, meta_data, wassr, delta_ppm, ppm_range))]
    
    @staticmethod
    def action(image, meta_data, wassr, delta_ppm=0.001, ppm_range=None):
        # Get the frequency information from the meta-data
        ppm = utils.get_ppm(meta_data)
        
        # Load the image data
        image = nibabel.load(image)
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

        nibabel.save(nibabel.Nifti1Image(B0_ppm, image.affine), wassr)
