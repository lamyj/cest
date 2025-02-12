import numpy

gamma = 267522187.44
gamma_bar = gamma/(2*numpy.pi)

# NOTE: B0 is never mentioned
B0 = 3 # T
w0_Hz = gamma_bar * B0
w0 = gamma * B0 # rad/s

Delta_omega_Hz = 20e3
Delta_omega_ppm = Delta_omega_Hz/w0_Hz*1e6

rf_offsets_Hz = numpy.linspace(-10000, +30000, 500)
rf_offsets_ppm = rf_offsets_Hz/w0_Hz*1e6

# NOTE: duration is never mentioned
duration = 10 # s

def M0(species_a, species_b):
    return numpy.array([0, 0, species_a.M0, 0, 0, species_b.M0])
