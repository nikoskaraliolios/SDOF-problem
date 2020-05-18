'''
Calculates the matrix with the inner products of BBPs

i,j and p index the BPs
h is the size of the timestep
c is the damping coefficitent of the SDOF system
k is the stiffness of the SDOF system
'''


import numpy as np
import scipy.integrate as integrate

import math

from .BernsteinPols import BP, dBP

# calculate inner products of BBPs. returns float


def integral_BP(i: int, j: int, p: int, h: float, c: float):
	return integrate.quad(lambda t: math.exp(c*t) * BP(t, i, p, h) * BP(t, j, p, h), 0 , h)[0]

# calculate inner products of derivatives of BBPs. returns float

def integral_dBP(i: int, j: int, p: int, h: float, c: float):
	return integrate.quad(lambda t: math.exp(c*t) * dBP(t, i, p, h) * dBP(t, j, p, h), 0 , h)[0]


# calculate the relevant submatrices of the matrix of inner products. returns p by p np.
# array of floats

def construct_B(p: int, h: float, c: float, k: float):
	B = np.zeros((p,p), np.float64)
	for i in range(p):
		for j in range(i,p):
			B[i,j] = k * integral_BP(i+1, j+1, p, h, c) - integral_dBP(i+1, j+1, p, h, c)
			B[j,i] = B[i,j]
	return B[1:-1,:2], np.linalg.inv(B[1:-1,2:])