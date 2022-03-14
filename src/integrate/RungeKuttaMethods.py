"""Runge-Kutta methods
"""

import numpy as np


class RungeKutta(object):
    """Base class for Runge-Kutta methods.
    """
    nstage: int = NotImplemented
    butcherTable: np.ndarray = NotImplemented

    def __init__(self):
        self.beta = self.butcherTable[:-1, 1:]
        self.w = self.butcherTable[-1, 1:]
        self.alpha = self.butcherTable[:-1,0]

    def step(self):
        """Perform one integration step.
        """
        raise NotImplementedError

    def nstep(self, fun, u, t0=0, dt=1., n=1):
        """Perform n times integration steps.

        Args:
            fun (callable): Right-hand side of the equation. du/dt = fun(t, u).
            u (np.array): Initial state.
            t0 (float): Initial time.
            dt (float): Timestep size.
            n (int, optional): Number of steps. Defaults to 1.

        Returns:
            np.array: u(t0 + dt * n)
        """
        self.t = t0
        for i in range(1, n+1):
            u = self.step(fun, u, t0, dt)
            self.t = t0 + dt*i
        return u


class RungeKutta2(RungeKutta):
    """Base class for Runge-Kutta methods of order 2.
    """
    def step(self, fun, u, t0=0, dt=1.):
        """Perform one integration step.

        Args:
            fun (callable): Right-hand side of the system.
            u (np.array): Initial state.
            t0 (float): Initial time.
            dt (float): Timestep size.

        Returns:
            np.array: u(t + dt)
        """
        k0 = dt * fun(t0, u)
        k1 = dt * fun(t0 + dt*self.alpha[1], u + self.beta[1,0]*k0)
        return u + self.w[0]*k0 + self.w[1]*k1


class RungeKutta3(RungeKutta):
    """Base class for Runge-Kutta methods of order 3.
    """
    def step(self, fun, u, t0=0, dt=1.):
        """Perform one integration step.

        Args:
            fun (callable): Right-hand side of the system.
            u (np.array): Initial state.
            t0 (float): Initial time.
            dt (float): Timestep size.

        Returns:
            np.array: u(t + dt)
        """      
        k0 = dt * fun(t0, u)
        k1 = dt * fun(t0 + dt*self.alpha[1], u + self.beta[1,0]*k0)
        k2 = dt * fun(t0 + dt*self.alpha[2], u + self.beta[2,0]*k0 + self.beta[2,1]*k1)
        return u + self.w[0]*k0 + self.w[1]*k1 + self.w[2]*k2


class RK2(RungeKutta2):
    """Runge-Kutta methods of order 2.

    The Butcher Table is from Table I Case no.2 [1]_.

    References:
    .. [1] J.H Williamson, "Low-storage Runge-Kutta schemes", 1980.
    """
    butcherTable = np.array([[0, 0, 0],
                             [2/3, 2/3, 0],
                             [0, 1/4, 3/4]])


class RK3(RungeKutta3):
    """Classic Runge-Kutta methods of order 3.

    The Butcher Table is from Table I Case no.11 [1]_.

    References:
    .. [1] J.H Williamson, "Low-storage Runge-Kutta schemes", 1980.
    """
    butcherTable = np.array([[0, 0, 0, 0],
                             [1/2, 1/2, 0, 0],
                             [1, -1, 2, 0],
                             [0, 1/6, 2/3, 1/6]])
