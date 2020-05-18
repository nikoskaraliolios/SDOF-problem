
'''
Calclulates numerically the indefinite integral of a function
given its sampling at a list of equally spaced points.

Takes the list of values of values and the sampling rate
(the distance between the points where the function is calculated)

Returns the values of the indefinite integral of the function
calculated at the same sequence of points by application of the
trapezoidal rule.
The value of the integral at the first point is 0.
'''


import numpy as np
from itertools import accumulate
import operator

def indefinite_integral(a: np.array((None,1), np.float), sampling_rate: float):
    trapezium = np.zeros_like(a)
    trapezium[1:] = (a[:-1]+a[1:]) * sampling_rate/2
    return np.asarray(list(accumulate(trapezium, operator.add)))
