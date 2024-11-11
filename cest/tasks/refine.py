import nibabel
import numpy
import scipy
import spire

from . import utils

class Refine(spire.TaskFactory):
    """Interpolate a Z-spectrum.
    
    Parameters
    ----------
    z_spectrum : path_like
        Path to the source Z-spectrum image
    meta_data : path_like
        Path to the meta-data related to the source image
    ppms : array_like
        Target frequencies
    refined : path_like
        Path to the target interpolated Z-spectrum
    """
    
    def __init__(self, z_spectrum, meta_data, ppms, refined):
        spire.TaskFactory.__init__(self, str(refined))
        self.file_dep = [z_spectrum, meta_data]
        self.targets = [refined]
        self.actions = [
            (__class__.action, (z_spectrum, meta_data, ppms, refined))]
    
    @staticmethod
    def action(z_spectrum, meta_data, ppms, refined):
        source_ppms = utils.get_ppm(meta_data)
        
        z_spectrum_image = nibabel.load(z_spectrum)
        z_spectrum = numpy.array(z_spectrum_image.dataobj)
        
        interpolator = scipy.interpolate.Akima1DInterpolator(
            source_ppms, z_spectrum, axis=3)
        refined_data = interpolator(ppms)
        
        nibabel.save(
            nibabel.Nifti1Image(refined_data, z_spectrum_image.affine), refined)

