"""Solvers of Advection Equation
"""

import numpy as np
from scipy.interpolate import CubicSpline

from integrate.lowStorageRungeKuttaMethods import LSRK3
from derivative.finiteDifferenceMethods import centralFDM, upwindFDM
from derivative.spectralMethods import spectralMethod1


class advection(object):
    def __init__(self, x, c=1.):
        self.x = x
        self.c = c
        self.nx = len(self.x)
        self.dx = np.diff(self.x)[0]
        self.integrate = LSRK3()

    def du_dx():
        raise NotImplementedError

    def rhs(self, t, u):
        return - self.c * self.du_dx(u)

    def solve(self, u0, t0, dt, ts):
        return self.integrate.nstep(fun=self.rhs, u=u0, t0=t0, dt=dt, n=ts)


class advectionCentralFDM(advection):
    def __init__(self, x, c=1., order=2):
        super().__init__(x=x, c=c)
        self.order = order
        self.derivative = centralFDM(order=self.order)
    
    def du_dx(self, u):
        return self.derivative(u, h=self.dx)


class advectionUpwindFDM(advection):
    def __init__(self, x, c=1., order=1):
        super().__init__(x=x, c=c)
        self.order = order
        self.derivative = upwindFDM(order=self.order)

    def du_dx(self, u):
        return self.derivative(u, h=self.dx, c=self.c)


class advectionSpectral(advection):
    def __init__(self, x, c=1.):
        super().__init__(x=x, c=c)
        spectralMethod = spectralMethod1(n=self.nx, d=self.dx)
        self.du_dx = spectralMethod.diff_phys


class advectionSciPyCubicSpline(advection):
    def __init__(self, x, c=1.):
        super().__init__(x=x, c=c)
        self.x_ = np.append(self.x, self.x[-1]+self.dx)

    def du_dx(self, u):
        u_ = np.append(u, u[0])
        cs = CubicSpline(self.x_, u_, bc_type='periodic')
        return cs(self.x_, nu=1, extrapolate='periodic')[:-1]


class advectionNumPyGradient(advection):
    def __init__(self, x, c=1.):
        super().__init__(x=x, c=c)

    def du_dx(self, u):
        u_ = np.hstack((u[-1], u, u[0]))
        u_ = np.gradient(u_, self.dx)
        return u_[1:-1]
