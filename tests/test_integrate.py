"""Test derivative methods
"""

import numpy as np

from cfdpy.derivative.spectralMethods import spectralMethod1
from cfdpy.integrate.RungeKuttaMethods import RK2, RK3
from cfdpy.integrate.lowStorageRungeKuttaMethods import LSRK2, LSRK3


def target_u(x):
    return np.sin(x) + np.sin(4*x) + 1.


def exact_u(fun, x, t, c=1.):
    return fun(x-c*t)


class advection(object):
    def __init__(self, x, c=1., integrate=None):
        self.x = x
        self.c = c
        self.nx = len(self.x)
        self.dx = np.diff(self.x)[0]

        spectralMethod = spectralMethod1(n=self.nx, d=self.dx)
        self.du_dx = spectralMethod.diff_phys

        self.integrate = integrate()

    def rhs(self, t, u):
        return - self.c * self.du_dx(u)

    def solve(self, u0, t0, dt, ts):
        return self.integrate.nstep(fun=self.rhs, u=u0, t0=t0, dt=dt, n=ts)


def test_integrator():
    c = 1.

    nx = 64
    x = np.linspace(0, 2*np.pi, nx, endpoint=False)

    t0 = 0. # Initial time.
    tf = np.abs(2*np.pi / c)
    nt = nx * 5
    dt = (tf - t0) / nt

    u0 = target_u(x)
    uf = exact_u(fun=target_u, x=x, t=tf, c=c)

    solver_rk2 = advection(x=x, c=c, integrate=RK2)
    solver_rk3 = advection(x=x, c=c, integrate=RK3)
    solver_lsrk2 = advection(x=x, c=c, integrate=LSRK2)
    solver_lsrk3 = advection(x=x, c=c, integrate=LSRK3)

    np.allclose(a=solver_rk2.solve(u0=u0, t0=t0, dt=dt, ts=nt), b=uf)
    np.allclose(a=solver_rk3.solve(u0=u0, t0=t0, dt=dt, ts=nt), b=uf)
    np.allclose(a=solver_lsrk2.solve(u0=u0, t0=t0, dt=dt, ts=nt), b=uf)
    np.allclose(a=solver_lsrk3.solve(u0=u0, t0=t0, dt=dt, ts=nt), b=uf)
