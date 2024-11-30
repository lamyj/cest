import numpy
import scipy

def two_pools(species_a, species_b, w, w1, duration, B0, Delta_B0=0):
    """ 
    Parameters
    ----------
    
    w: float
        frequency offset of the saturation RF pulse in Hz
    w1: float
        frequency of the B1 field of the saturation RF pulse in Hz
    duration: float
        duration in s
    B0: float
        nominal B0 field in Hz
    Delta_B0: float, optional
        B0 offset in PPM
    """
    R1a, R2a = 1/species_a.T1, 1/species_a.T2
    R1b, R2b = 1/species_b.T1, 1/species_b.T2
    wa, wb = species_a.w*1e-6 * B0, species_b.w*1e-6 * B0
    M0a, M0b = species_a.M0, species_b.M0
    
    fb = species_b.M0/species_a.M0
    
    ka = species_b.k
    
    w0 = Delta_B0*1e-6 * B0
    
    pi = numpy.pi
    
    omega_a = 2*pi*(w-wa+w0)
    omega_b = 2*pi*(w-wb+w0)
    
    fbka = fb*ka
    
    A = numpy.array([
        [-R2a-fbka,  -omega_a,         0,      ka,        0,       0,       0],
        [  omega_a, -R2a-fbka,       -w1,       0,       ka,       0,       0],
        [        0,        w1, -R1a-fbka,       0,        0,      ka, R1a*M0a],
        [     fbka,         0,         0, -R2b-ka, -omega_b,       0,       0],
        [        0,      fbka,         0, omega_b,  -R2b-ka,     -w1,       0],
        [        0,         0,      fbka,       0,       w1, -R1b-ka, R1b*M0b],
        [        0,         0,         0,       0,        0,       0,       0]])
      
    return scipy.linalg.expm(A*duration)
