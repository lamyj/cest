import base64
import json

import dicomifier
import numpy

def get_ppm(path):
    """
    Return the frequency offsets as PPM from the meta-data.
    
    Parameters
    ----------
    path : path_like
        Path to the meta-data
    """
    with open(path) as fd:
        meta_data = json.load(fd)
    
    encapsulated = base64.b64decode(meta_data["EncapsulatedDocument"][0])
    if encapsulated[-1] == 0:
        encapsulated = encapsulated[:-1]
    data = json.loads(encapsulated)
    
    data_set = dicomifier.bruker.Dataset()
    data_set.loads(data["acqp"])
    data_set.loads(data["method"])
    
    ppm = None
    
    omega_0_MHz = data_set["BF1"].value[0]
    if "PVM_MagTransFL" in data_set:
        delta_omega_Hz = data_set["PVM_MagTransFL"].value
        ppm = numpy.divide(delta_omega_Hz, omega_0_MHz)
    # TODO: Otherwise use CESTFreqList or PVM_MagTransOffset
    else:
        raise Exception("Cannot find PPM information")
    
    return ppm
