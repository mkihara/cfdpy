# CFDPy

Computational Fluid Dynamics for Python.

## Installation

```sh
pip install cfdpy
```

## Features

There are already some usefull packages for scientific computing with Python, such as

* [NumPy](https://numpy.org/)
* [SciPy](https://scipy.org/)

Based on these packages, CFDPy has implemented some additional features for CFD.
CFDPy is different from the previous packages in the following ways:

1. CFDPy's ODE solvers accepts multidimensional arrays and adopts faster (low-accuracy) methods.
2. CFDPy's derivative algorithm supports higher order differentials and accuracy.

CFDPy focuses on simplicity and is not optimized for parallel computing. It is useful for learning CFD, but not for fast computation.
If you have faster computation needs, we recommend exploring Fortran, C++ or Julia codes.

## Current functionality

* Finite difference methods of arbitrary order differentials and accuracy.
* De-aliasing spectral methods.
* Runge-Kutta methods of order 2 and 3.
* Low-storage Runge-Kutta methods of order 2 and 3.
