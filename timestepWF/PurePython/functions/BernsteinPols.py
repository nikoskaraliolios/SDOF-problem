
'''
Define Bernstein Polynomials

i and p index the BPs
t is the time variable
h is the size of the timestep
'''


import scipy.special
import math




def BP(t: float, i: int, p: int, h: float):
	if i < 1 or i > p:
		return 0
	elif t >= 0 and t <= h:
		norm_coef = scipy.special.binom(p-1,i-1)
		t = t/h
		return norm_coef * t**(i-1) *  (1-t)**(p-i)
	else:
		return 0



def dBP(t: float, i: int, p: int, h: float):
	if t >= 0 and t <= h:
		return (p-1)*(BP(t, i-1, p-1, h) - BP(t, i, p -1, h) ) / h
	else:
		return 0