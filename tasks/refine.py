import nibabel
import numpy
import scipy
import spire

from . import utils

class Refine(spire.TaskFactory):
    def __init__(self, z_spectrum, meta_data, frequencies, refined):
        spire.TaskFactory.__init__(self, str(refined))
        self.file_dep = [z_spectrum, meta_data]
        self.targets = [refined]
        self.actions = [
            (__class__.action, (z_spectrum, meta_data, frequencies, refined))]
    
    @staticmethod
    def action(z_spectrum_path, meta_data_path, frequencies, refined_path):
        ppm = utils.get_ppm(meta_data_path)
        
        z_spectrum_image = nibabel.load(z_spectrum_path)
        z_spectrum = numpy.array(z_spectrum_image.dataobj)
        
        interpolator = scipy.interpolate.Akima1DInterpolator(
            ppm, z_spectrum, axis=3)
        refined = interpolator(frequencies)
        
        nibabel.save(
            nibabel.Nifti1Image(refined, z_spectrum_image.affine), refined_path)
