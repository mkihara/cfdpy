"""Test derivative methods
"""

import numpy as np

from cfdpy.derivative.finiteDifferenceMethods import centralFDM, upwindFDM
from cfdpy.derivative.spectralMethods import spectralMethod1


def target_u(x):
    return np.sin(x) + np.sin(4*x) + 1.


def exact_du(x):
    return np.cos(x) + 4*np.cos(4*x)


def exact_ddu(x):
    return - np.sin(x) - 16*np.sin(4*x)


def test_finiteDifferenceMethods():
    nx = 64
    x = np.linspace(0, 2*np.pi, nx, endpoint=False)
    dx = 2*np.pi / nx
    u = target_u(x)
    du = exact_du(x)
    ddu = exact_ddu(x)

    cfdm2 = centralFDM(order=2, highestDerivative=2)
    cfdm4 = centralFDM(order=4, highestDerivative=2)
    cfdm6 = centralFDM(order=6, highestDerivative=2)
    np.allclose(a=cfdm2(f=u, derivative=1, h=dx), b=du)
    np.allclose(a=cfdm4(f=u, derivative=1, h=dx), b=du)
    np.allclose(a=cfdm6(f=u, derivative=1, h=dx), b=du)
    np.allclose(a=cfdm2(f=u, derivative=2, h=dx), b=ddu)
    np.allclose(a=cfdm4(f=u, derivative=2, h=dx), b=ddu)
    np.allclose(a=cfdm6(f=u, derivative=2, h=dx), b=ddu)

    ufdm1 = upwindFDM(order=1, highestDerivative=2)
    ufdm3 = upwindFDM(order=3, highestDerivative=2)
    ufdm5 = upwindFDM(order=5, highestDerivative=2)
    np.allclose(a=ufdm1(f=u, derivative=1, h=dx), b=du)
    np.allclose(a=ufdm3(f=u, derivative=1, h=dx), b=du)
    np.allclose(a=ufdm5(f=u, derivative=1, h=dx), b=du)
    np.allclose(a=ufdm1(f=u, derivative=2, h=dx), b=ddu)
    np.allclose(a=ufdm3(f=u, derivative=2, h=dx), b=ddu)
    np.allclose(a=ufdm5(f=u, derivative=2, h=dx), b=ddu)


def test_spectralMethods():
    nx = 64
    x = np.linspace(0, 2*np.pi, nx, endpoint=False)
    dx = 2*np.pi / nx
    u = target_u(x)
    du = exact_du(x)

    spectral1 = spectralMethod1(n=nx, d=dx)
    np.allclose(a=spectral1.diff_phys(u), b=du)
