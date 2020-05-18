<p align="center"><img src="/tex/ceec6b2a02661b264e83a9dce4de6d39.svg?invert_in_darkmode&sanitize=true" align=middle width=402.4065111pt height=17.5342167pt/></p>

This repository contains an implementation of the weak solution of the SDOF (Single Degree of Freedom) problem in Python.

The SDOF problem consists in solving the problem
<p align="center"><img src="/tex/8fdf16d6db753f313f17a88e2b54421b.svg?invert_in_darkmode&sanitize=true" align=middle width=130.724616pt height=14.611878599999999pt/></p>
with arbitrary initial conditions
<p align="center"><img src="/tex/0985d062e65bd9b9ab3c8dc4b304a190.svg?invert_in_darkmode&sanitize=true" align=middle width=72.80807325pt height=14.52054615pt/></p>
and with an arbitrary excitation function
<p align="center"><img src="/tex/56dfde2ff4e655b1b00c0d9cc46e2a46.svg?invert_in_darkmode&sanitize=true" align=middle width=97.5054432pt height=17.5981443pt/></p>
acting on the system for a finite time.


This is relevant to Civil Engineering problems, as it modelizes the response of simple
structures under a dynamic load (e.g. an earthquake or a train crossing a bridge).
More relevant, of course, is the MDOF (Multi-Degree of Freedom) problem, which models
more compex structures.


The theoretical background generating a family of algorithms for both the SDOF and
MDOF problems will be provided in a forthcoming article in collaboration with prof.
Dimitrios L. Karabalis. It will be followed up by an implementation of the general
framework in a numerical method based on a time-step polynomial approximation scheme.

Numerical experiments carried out using Matlab (written by the author) and Fortran
(coded by prof. Stefanos Tsinopoulos)
already show that the algorithm proposed in the forthcoming article improves the state
of the art by several orders of magnitude, without significantly deteriorating the
computational cost.

The subfolder "WF in each timestep" contains a pure Python implementation of that method,
accompanied by Jupyter notebooks for performing some basic tests.

The authors, N. Karaliolios and D.L. Karabalis, will come back to the problem and
provide a numerical method for the solution of the MDOF problem.

TODOs: The repository will be updated with a Cythonized and a parallelized version of the
same algorithm, as well as with a second, more efficient numerical method.
