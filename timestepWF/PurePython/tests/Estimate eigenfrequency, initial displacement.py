
'''
Estimate the eigenfrequency the SDOF problem through changes
of direction of motion when x0 =1, and dx0 and f are 0

Displays a plot of the displacement and of the speed, and
prints out the error in the estimation of the half eigenperiods
contained in the interval of approximation
'''

import sys
sys.path.append('../')


import numpy as np
import math

import matplotlib.pyplot as plt

from functions.BernsteinPols import BP, dBP
from functions.construct_B import construct_B
from functions.iteration import SDOF_integration
from functions.make_solution import get_displacement, get_speed, get_acceleration, get_solution


# Define parameters for the SDOF system



m   = 1.5               # the mass, >0
c   = 0.1               # the damping coefficient, >= 0
k   = 2                 # the stiffness coefficient, >0
x0  = 1                 # the initial displacement
dx0 = 0                 # the initial speed
def force(t):           # the excitation function = 0 for this test
	return 0
    # return np.sin(t)


# Define parameters for the approximation



p   = 5                    # p-1 is the degree of polynomial approximation, p>= 3
h   = 0.1                  # the timestep, >0
l   = 100                  # the number of iterations, >0
sampling_rate = 0.00001    # the spacing of points to sample for estimating the eigenfrequency


# Initialize and run the algorithm



c = c/m
k = k/m
def f(t):
    return force(t)/m

sdof_int = SDOF_integration(p, c, k, x0, dx0, f, h, l)

BP_coefs = sdof_int.integrator()
t_series, approx_disp, approx_speed, _ = get_solution(BP_coefs, c, k, f, h, sampling_rate = sampling_rate)
approx_speed = approx_speed[:t_series.shape[0]]


# Calculate and plot the approximate speed




fig = plt.figure()
# ax = plt.axes()
ax = fig.gca()
ax.set_xticks(np.arange(-0, l*h, 1))
ax.set_yticks(np.arange(-1., 1.1, 0.5))
plt.title("Approximate solution")
ax.plot(t_series, approx_speed, 'red', label = 'Speed');
ax.plot(t_series, approx_disp, 'blue', label = 'Displacement');
ax.legend();
plt.grid()
plt.show()


# Find changes of sign in speed



change_direction = list(map(lambda x:x[0]*x[1] <= 0,list(zip(approx_speed[:-1],approx_speed[1:]))))
idx = [sampling_rate*i for i,x in enumerate(change_direction) if x]
xi = c / (2*np.sqrt(k))
true_period = np.pi / (np.sqrt(k) * np.sqrt(1 - xi **2))
error = (idx[1:] - true_period*np.arange(1,len(idx)))/ true_period

for i in range(len(error)):
	if i == 0:
		print("The relative error in the estimation of the {i}st half-eigenperiod is {e:.2e}".format(i = i+1, e = error[i]))

	elif i == 1:
		print("The relative error in the estimation of the {i}nd half-eigenperiod is {e:.2e}".format(i = i+1, e = error[i]))

	elif i == 2:
		print("The relative error in the estimation of the {i}rd half-eigenperiod is {e:.2e}".format(i = i+1, e = error[i]))

	else:
		print("The relative error in the estimation of the {i}th half-eigenperiod is {e:.2e}".format(i = i+1, e = error[i]))

