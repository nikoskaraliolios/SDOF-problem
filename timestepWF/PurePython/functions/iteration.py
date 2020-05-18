'''
The class SDOF_integration contains the information of:
i) the SDOF system
ii) the paramters of approximation


Respectively:
i) c, k, the exictation function f and
the initial conditions x0 and dx0; NB MASS IS NORMALIZED TO 1
ii) the parameters of the approximation scheme h and l


The integrator function takes SDOF_integration object as its only
argument, and returns the p by l np.array with the coefficients of
the BP that approximate the solution to the problem corresponding
to the SDOF_integration object.
'''


import numpy as np
import math

from .construct_B import construct_B
from .construct_rhs import construct_rhs


class SDOF_integration():
	def __init__(self, p: int, c: float, k: float, x0: float, dx0: float, f, h: float, l: int):
		self.p   = p							# the degree of polynomial approximation, >= 3
		self.c   = c							# the damping coefficient
		self.k   = k							# the stiffness coeffient
		self.x0  = x0							# the excitation function
		self.dx0 = dx0							# the stiffness coeffient
		self.f   = f							# the excitation function
		self.h   = h							# the timestep
		self.l   = l							# the number of timesteps


		self.cols, self.coeff_mat = construct_B(p, h, c, k)
		self.F = construct_rhs(p, f, h, c, k, l)
		self.s = h/(p-1)


	'''
	solves the linear system corresponding to one timestep
	cols and coef_mat are the output of construct_B
	F_step is one element of the output of construct_rhs along axis = 0
	first_two_coefs are the first two BP coefficients for the timestep,
	which depend directly on the initial displacement and speed.
	
	returns an array of BP coefficients for the corresponding timestep.
	'''
	@staticmethod
	def _step_solver(cols: np.array, coef_mat: np.array, F_step: np.array, first_two_coefs: np.array):
		F_step -= np.dot(cols, first_two_coefs)
		return np.dot(coef_mat, F_step)
	
	def step_solver(self, F_step: np.array, first_two_coefs: np.array):
		return self._step_solver(self.cols, self.coeff_mat, F_step, first_two_coefs)

	'''
	Takes an instance of SDOF_integration and returns the BP coefficients that
	construct the approximate solution
	'''

	def integrator(self):
		BP_coefs = np.empty((self.p,self.l+1), np.float)
		BP_coefs[0,0]  = self.x0
		BP_coefs[1,0]  = self.s * self.dx0 + self.x0
		
		for q in range(self.l):
			BP_coefs[2:,q] = self.step_solver(self.F[:,q],BP_coefs[:2,q])
			BP_coefs[0,q+1] = BP_coefs[-1,q]
			BP_coefs[1,q+1] = 2*BP_coefs[-1,q] - BP_coefs[-2,q]
		return BP_coefs[:,:-1]