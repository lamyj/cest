"""
Library of RF pulses. Each pulse function has the same first parameter, i.e. the
number of steps of the pulse support. The integral of each pulse is normalized
"""

import numpy

def block(steps):
    y = numpy.ones(steps)
    return y/y.sum()
    
def gaussian(steps, sigma=1, x_max=3.5):
    x = numpy.linspace(-x_max, x_max, steps)
    y = numpy.exp(-(x**2 / (2*sigma**2))) / (sigma * numpy.sqrt(2*numpy.pi))
    return y/y.sum()

def sinc(steps, side_lobes=2):
    x = numpy.linspace(-side_lobes-1, +side_lobes+1, steps)
    y = numpy.sinc(x)
    return y/y.sum()

def sech(steps, x_max=10):
    x = numpy.linspace(-1, +1, steps)
    y = 1/numpy.cosh(x_max * x)
    return y/y.sum()

def train(pulse, count, gap_steps):
    y = numpy.empty(count*len(pulse)+(count-1)*gap_steps)
    begin = 0
    for i in range(count-1):
        y[begin:begin+len(pulse)] = pulse
        begin += len(pulse)
        
        y[begin:begin+gap_steps] = 0
        begin += gap_steps
    
    y[-len(pulse):] = pulse
    
    return y/y.sum()
