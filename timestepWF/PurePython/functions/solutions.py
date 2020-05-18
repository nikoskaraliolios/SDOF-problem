'''
Provides the exact solutions to some SDOF problems, in order
to run tests and estimate the efficiency of the method

All functions return displacement, speed and acceleration in
this order
'''

import numpy as np


# arbitrary parameters and initial conditions, 0 excitation function

def hom_sol(t: float, c: float, k: float, x0: float, dx0: float):
	
	w_n = np.sqrt(k)
	
	xi  = c/(2 * w_n)

	w_d = w_n * np.sqrt(1 - xi **2)

	exact_disp = np.exp(- xi * w_n * t)*(x0 * np.cos(w_d *t) +
						(dx0 + xi * w_n * x0)/w_d * np.sin(w_d *t))
		
	aux_term = - w_d *x0 * np.sin(w_d *t) + (dx0 + xi * w_n * x0) * np.cos(w_d *t)
	d_aux_term = - w_d **2 *x0 * np.cos(w_d *t) - (dx0 + xi * w_n * x0) * w_d * np.sin(w_d *t)

	exact_speed = - xi * w_n * exact_disp + np.exp(- xi * w_n * t)* aux_term

	exact_acc =  - xi * w_n * exact_speed - xi * w_n * np.exp(- xi * w_n * t)* aux_term + np.exp(- xi * w_n * t)* d_aux_term

	return exact_disp, exact_speed, exact_acc

# m = k = 1, c = 0, sin(t) excitation function, 0 initial conditions

def sin_exc_fun(t: float):
	return (np.sin(t) - t* np.cos(t))/2, t * np.sin(t)/2, (np.sin(t)+ t* np.cos(t))/2


# m = k = 1, c = 0, sin(2t) excitation function, 0 initial conditions

def sin2_exc_fun(t: float):
	return -2*np.sin(t) *(np.cos(t)-1)/3, 2*(np.cos(t) - np.cos(2*t))/3, -2*(np.sin(t) - 2* np.sin(2*t))/3


