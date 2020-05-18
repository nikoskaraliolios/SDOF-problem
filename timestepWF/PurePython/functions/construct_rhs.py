
'''
Calculates the inner products of
the excitation function f with BPs

i and p index the BPs
f is the excitation function
h is the size of the timestep
c is the damping coefficitent of the SDOF system
q is the index of the time step, in [0,l-1]
'''


import numpy as np
import scipy.integrate as integrate

import math

from .BernsteinPols import BP, dBP

# calculates the the inner product of f with the i-th BP translated by q timesteps
# returns float

def integral_F(i: int, f, p: int, h: float, c: float, q: int):
	return integrate.quad(lambda t: math.exp(c*t) * BP(t , i, p, h) * f(t +q*h), 0 , h)[0]


# calculate the vector of the inner products of f with BPs for each timestep
# returns p-2 by 1 np array of floats

def construct_F_step(p: int, f, h: float, c: float,k: float, q: int):
	F_step = np.zeros((p-2), np.float64)
	for i in range(2,p):
		F_step[i-2] = integral_F(i, f, p, h, c, q)
	
	return F_step

# calculate the vector for all timesteps
# returns p-2 by l np array of floats

def construct_rhs(p: int, f, h: float, c: float,k: float, l: int):
	F = np.zeros((p-2, l), np.float64)
	for i in range(l):
		F[:,i] = construct_F_step(p, f, h, c, k, i)
	
	return F