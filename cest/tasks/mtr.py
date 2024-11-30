import pathlib

import nibabel
import numpy
import spire

from . import utils

class MTR(spire.TaskFactory):
    """Compute an MTR map.
        
    Parameters
    ----------
    
    z_spectrum : path-like or image or array
        Source Z-spectrum data. It must be path-like if used as a task,
        otherwise it may be a nibabel image or a numpy array.
    ppms : path-like or array
        PPM values of the z-spectrum. If path-like, the values will be loaded
        from the path. The PPM values must be symmetric with respect to 0.
    mtr : path-like, optional
        Path to the target MTR image. It is optional only when used as a
        function.
    normalization : string, optional
        Normalization method, must be one of *asym*, *normref* (default 
        value), *pcm*, *rex*
    
    Returns
    -------
    
    array or image
        MTR data, only applicable when used as a function.
    
    References
    ----------
    
    *Inverse Z-spectrum analysis for spillover-, MT-, and T1-corrected
    steady-state pulsed CEST-MRI – application to pH-weighted MRI of acute
    stroke*. Zaiss et al. NMR in Biomedicine 27(3), 2014.
    `doi:10.1002/nbm.3054 <https://doi.org/10.1002/nbm.3054>`_.
    """
    def __init__(self, z_spectrum, ppms, mtr, normalization="normref"):
        spire.TaskFactory.__init__(self, str(mtr))
        self.file_dep = [
            z_spectrum,
            *([ppms] if isinstance(ppms, (str, pathlib.Path)) else [])]
        self.targets = [mtr]
        self.actions = [
            (__class__.action, (z_spectrum, ppms, mtr, normalization))]
    
    @staticmethod
    def action(z_spectrum, ppms, mtr=None, normalization="normref"):
        if isinstance(ppms, (str, pathlib.Path)):
            ppms = utils.get_ppm(ppms)
        
        # WARNING: must be symmetrical, as stated in docstring
        positive = ppms >=0
        negative = ppms <= 0
        
        if isinstance(z_spectrum, (str, pathlib.Path)):
            z_spectrum = nibabel.load(z_spectrum)
        if isinstance(z_spectrum, nibabel.Nifti1Image):
            image = z_spectrum
            z_spectrum = numpy.asarray(z_spectrum.dataobj)
            affine = image.affine
        else:
            image = None
            affine = numpy.identity(4)
        
        z_label = z_spectrum[..., positive]
        z_reference = z_spectrum[..., negative][::-1]
        
        if normalization.lower() == "asym": # Eq. 7
            mtr_data = z_reference - z_label
        elif normalization.lower() == "normref": # Eq. 8
            mtr_data = (z_reference - z_label) / z_reference
        elif normalization.lower() == "pcm": # Eq. 9
            mtr_data = (
                (z_reference - z_label)
                / (z_reference - z_label + z_label*z_reference))
        elif normalization.lower() == "rex": # Eq. 10
            mtr_data = 1/z_label - 1/z_reference
        else:
            raise Exception(f"Unknown normalization: {normalization}")
        
        if mtr is None:
            return (
                nibabel.Nifti1Image(mtr_data, affine) if image is not None
                else mtr_data)
        
        nibabel.save(nibabel.Nifti1Image(mtr_data, affine), mtr)

