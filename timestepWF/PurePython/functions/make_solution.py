'''
All functions take as arguments
i) a p by l np.array of floats. They are the
coefficients of the BP as returned by the
integrator function of .make_solution
ii) the timestep h>0
iii) the sampling_rate >0. This is the spacing
between the points where the solutions are sampled
in order to plot the functions and perform the tests.

All functions return a tuple of np arrays. The first
is the time-series of sampling points and the rest
are the corresponding values of the (displacement, speed
and/or acceleration) of the solution at the points
of the time-series
'''




import numpy as np
import math

from .BernsteinPols import BP, dBP



def get_displacement(BP_coefs: np.array, h: float, sampling_rate = 0.001):
	
	p = BP_coefs.shape[0]
	l = BP_coefs.shape[1]

		
	sample_one_step = np.arange(0, h, sampling_rate)
	sample = np.arange(0,l * h, sampling_rate)

	
	BP_sample = np.zeros((p, sample_one_step.shape[0]))
	
	for i in range(p):
		BP_sample[i,:] = [BP(t,i+1,p,h) for t in sample_one_step]
	
	
	approx_disp =  np.zeros((sample_one_step.shape[0],l))
	
	for j in range(l):
		approx_disp[:,j] = np.dot(BP_coefs[:,j],BP_sample)
	
	approx_disp = approx_disp.flatten("F")
	approx_disp = approx_disp[:sample.shape[0]]
	
	print(approx_disp.shape)

	return sample, approx_disp


def get_speed(BP_coefs: np.array, h: float, sampling_rate = 0.001):
	
	p = BP_coefs.shape[0]
	l = BP_coefs.shape[1]

		
	sample_one_step = np.arange(0, h, sampling_rate)
	sample = np.arange(0,l * h, sampling_rate)

	
	dBP_sample = np.zeros((p, sample_one_step.shape[0]))
	
	for i in range(p):
		dBP_sample[i,:] = [dBP(t,i+1,p,h) for t in sample_one_step]
	
	
	approx_speed =  np.zeros((sample_one_step.shape[0],l))
	
	for j in range(l):
		approx_speed[:,j] = np.dot(BP_coefs[:,j],dBP_sample)
	
	approx_speed = approx_speed.flatten("F")
	approx_speed = approx_speed[:sample.shape[0]]
	
	return sample, approx_speed


def get_acceleration(BP_coefs: np.array, c, k, f, h: float, sampling_rate = 0.001):
	
	# p = BP_coefs.shape[0]	
	l = BP_coefs.shape[1]

	
	sample = np.arange(0,l * h, sampling_rate)
	f_sample = [f(t) for t in sample]

	_, approx_disp = get_displacement(BP_coefs, h)
	_, approx_speed = get_speed(BP_coefs, h)
	
	return sample, f_sample - c * approx_speed - k * approx_disp



def get_solution(BP_coefs: np.array, c, k, f, h: float, sampling_rate = 0.001):
	
	p = BP_coefs.shape[0]
	l = BP_coefs.shape[1]

		
	sample_one_step = np.arange(0, h, sampling_rate)
	sample = np.arange(0,l * h, sampling_rate)

	BP_sample = np.zeros((p, sample_one_step.shape[0]))
	dBP_sample = np.zeros((p, sample_one_step.shape[0]))
	f_sample = [f(t) for t in sample]
	
	for i in range(p):
		BP_sample[i,:]  = [ BP(t,i+1,p,h) for t in sample_one_step]
		dBP_sample[i,:] = [dBP(t,i+1,p,h) for t in sample_one_step]
	
	approx_disp =  np.zeros((sample_one_step.shape[0],l))
	approx_speed =  np.zeros((sample_one_step.shape[0],l))
	
	for j in range(l):
		approx_disp[:,j] = np.dot(BP_coefs[:,j],BP_sample)
		approx_speed[:,j] = np.dot(BP_coefs[:,j],dBP_sample)

	
	approx_disp = approx_disp.flatten("F")
	approx_disp = approx_disp[:sample.shape[0]]
	approx_speed = approx_speed.flatten("F")
	approx_speed = approx_speed[:sample.shape[0]]

	
	return sample, approx_disp, approx_speed, f_sample - c * approx_speed - k * approx_disp
