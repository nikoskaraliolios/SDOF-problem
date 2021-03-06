'''
Plot approximation errors of solutions to the homogeneous SDOF problem
with arbitrary parameters and initial conditions
'''

import sys
sys.path.append('../')

import numpy as np
import math

import matplotlib.pyplot as plt
import time

from functions import BernsteinPols
from functions.construct_B import construct_B
from functions.iteration import SDOF_integration
from functions.make_solution import get_displacement, get_speed, get_acceleration, get_solution


# Define parameters for the SDOF system



m   = 1.5                 # the mass, >0
c   = 0.1                 # the damping coefficient, >= 0
k   = 2                 # the stiffness coefficient, >0
x0  = 0.1                 # the initial displacement
dx0 = -0.1                 # the initial speed
def force(t):           # the excitation function
	return 0
    # return np.sin(t)


# Provide exact solution



from functions.solutions import hom_sol as exact

def sol(t):
    return exact(t, c, k, x0, dx0)

def exact_disp(t):
    return sol(t)[0]

def exact_speed(t):
    return sol(t)[1]

def exact_acc(t):
    return sol(t)[2]


# Define parameters for the approximation



p   = 5                 # p-1 is the degree of polynomial approximation, p>= 3
h   = 0.1               # the timestep, >0
l   = 100               # the number of iterations, >0


# Initialize and run the algorithm



t0 = time.time()
c = c/m
k = k/m
def f(t):
    return force(t)/m

sdof_int = SDOF_integration(p, c, k, x0, dx0, f, h, l)
t_init = time.time()
BP_coefs = sdof_int.integrator()
t_series, approx_disp, approx_speed, approx_acc = get_solution(BP_coefs, c, k, f, h)
t_solve = time.time() - t_init


# In[6]:


print("Initialization took {i:.2e} sec, and execution took {e:.2e} sec".format(i = t_init-t0, e = t_solve))


# Calculate and plot the error of approximation in displacement



error_disp = exact_disp(t_series) - approx_disp

fig = plt.figure()
ax = plt.axes()
plt.title("Error in displacement")
ax.plot(t_series, error_disp);
# plt.show()


# Calculate and plot the error of approximation in speed



error_speed = exact_speed(t_series) - approx_speed

fig = plt.figure()
ax = plt.axes()
plt.title("Error in speed")
ax.plot(t_series, error_speed);
# plt.show()


# Calculate and plot the error of approximation in acceleration



error_acc = exact_acc(t_series) - approx_acc

fig = plt.figure()
ax = plt.axes()
plt.title("Error in acceleration")
ax.plot(t_series, error_acc);
plt.show()

