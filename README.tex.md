\begin{center}
\begin{large}
\textbf{Numerical solution of the SDOF problem}
\end{large}
\end{center}

This repository contains an implementation of the weak solution of the SDOF (Single Degree of Freedom) problem in Python.

The SDOF problem consists in solving the problem
$$
m \ddot{x} + c \dot{x} + k x = f
$$
with arbitrary initial conditions
$$x_{0} , \dot{x}_{0} \in \mathbb{R}$$
and with an arbitrary excitation function
$$f : [0 , \bar{T}] \rightarrow \mathbb{R}$$
acting on the system for a finite time.


This is relevant to Civil Engineering problems, as it modelizes the response of simple
structures under a dynamic load (e.g. an earthquake or a train crossing a bridge).
More relevant, of course, is the MDOF (Multi-Degree of Freedom) problem, which models
more compex structures.


The theoretical background generating a family of algorithms for both the SDOF and
MDOF problems will be provided in a forthcoming article in collaboration with prof.
Dimitrios L. Karabalis. It will be followed up by an implementation of the general
framework in a numerical method based on a time-step polynomial approximation scheme.

Numerical experiments show that the algorithm proposed in the forthcoming article
improves the state of the art by several orders of magnitude, without deteriorating
the computational cost.

The subfolder "WF in each timestep" contains a pure Python implementation of that method,
accompanied by Jupyter notebooks for performing some basic tests.

The authors, N. Karaliolios and D.L. Karabalis, will come back to the problem and
provide a numerical method for the solution of the MDOF problem.

TODOs: The repository will be updated with a Cythonized and a parallelized version of the
same algorithm, as well as with a second, more efficient numerical method. Scripts and
notebooks performing some tests will also be added.
