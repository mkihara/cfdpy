"""Solvers of Burgers Equations
"""

import numpy as np

from integrate.lowStorageRungeKuttaMethods import LSRK3
from derivative.spectralMethods import spectralMethod1
from derivative.finiteDifferenceMethods import centralFDM


class burgers1(object):
    def __init__(self, x, nu=1.):
        self.x = x
        self.nu = nu
        self.nx = len(self.x)
        self.dx = np.diff(self.x)[0]
        self.integrate = LSRK3()
        self.spectral = spectralMethod1(n=self.nx, d=self.dx)

    def nonlinear(self, u):
        return self.spectral.multiply(u, self.spectral.diff(u))

    def viscosity(self, u):
        return self.nu * self.spectral.diff(u, order=2)

    def rhs(self, t, u):
        return self.viscosity(u) - self.nonlinear(u)

    def solve(self, u0, t0, dt, ts):
        u = self.spectral.fft(u0)
        u = self.integrate.nstep(fun=self.rhs, u=u, t0=t0, dt=dt, n=ts)
        return self.spectral.ifft(u).real


class burgers2(object):
    def __init__(self, dx, dy, nu=1.):
        self.dx = dx
        self.dy = dy
        self.nu = nu
        self.integrate = LSRK3()
        self.derivative = centralFDM(order=6, highestDerivative=2)

    def rhs(self, t, u):
        return ( - u[0] * self.derivative(u, h=self.dx, axis=-1)
                 - u[1] * self.derivative(u, h=self.dy, axis=-2)
                 + self.nu * (self.derivative(u, h=self.dx, axis=-1, derivative=2) 
                            + self.derivative(u, h=self.dy, axis=-2, derivative=2)) )

    def solve(self, u0, t0, dt, ts):
        return self.integrate.nstep(fun=self.rhs, u=u0, t0=t0, dt=dt, n=ts)
