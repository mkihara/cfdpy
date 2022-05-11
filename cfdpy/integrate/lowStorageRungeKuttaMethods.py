"""Low-storage Runge-Kutta methods
"""

import numpy as np

from .RungeKuttaMethods import RungeKutta


class lowStorageRungeKutta(RungeKutta):
    """Base class for low-storage Runge-Kutta methods.

    This uses the Williamson's 2N-storage scheme [1]_.

    References:
    .. [1] J.H Williamson, "Low-storage Runge-Kutta schemes", 1980.
    """
    def __init__(self):
        super().__init__()

        self.b = np.zeros(self.nstage)
        for j in range(self.nstage-1):
            self.b[j] = self.beta[j+1,j]
        self.b[self.nstage-1] = self.w[self.nstage-1]

        self.a = np.zeros(self.nstage)
        for j in range(1, self.nstage):
            if self.w[j] != 0:
                self.a[j] = (self.w[j-1] - self.b[j-1]) / self.w[j]
            else:
                self.a[j] = (self.beta[j+1,j-1] - self.alpha[j]) / self.b[j]

    def step(self, fun, u, t0=0, dt=1.):
        """Perform one integration step.

        Args:
            fun (callable): Right-hand side of the system.
            u (np.array): Initial state.
            t0 (float): Initial time.
            dt (float): Timestep size.

        Returns:
            np.array: u(t0 + dt)
        """
        ### j = 0
        du = dt * fun(t0+dt*self.alpha[0], u)
        u = u + self.b[0]*du
        ### j = 1 to n-1
        for j in range(1, self.nstage):
            du = self.a[j]*du + dt*fun(t0+dt*self.alpha[j], u)
            u += self.b[j]*du
        return u


class LSRK2(lowStorageRungeKutta):
    """Low-storage Runge-Kutta method of order 2.

    This uses the Williamson's 2N-storage scheme [1]_.
    The Butcher Table is from Table I Case no.2 [1]_.

    Attributes:
        nstage (int): Number of stage, which is equivalent to the order of accuracy.
        butcherTable (numpy.ndarray): [[alpha, beta], [0, w]]

    Methods:
        step(fun, u, t0, dt) : Perform one integration step.
        nstep(fun, u, t0, dt, nt) : Perform n times integration steps.

    References:
    .. [1] J.H Williamson, "Low-storage Runge-Kutta schemes", 1980.
    """
    nstage = 2
    butcherTable = np.array([[0, 0, 0],
                             [2/3, 2/3, 0],
                             [0, 1/4, 3/4]])


class LSRK3(lowStorageRungeKutta):
    """Low-storage Runge-Kutta method of order 3.

    This uses the Williamson's 2N-storage scheme [1]_.
    The Butcher Table is from Table I Case no.7 [1]_.

    Attributes:
        nstage (int): Number of stage, which is equivalent to the order of accuracy.
        butcherTable (numpy.ndarray): [[alpha, beta], [0, w]]

    Methods:
        step(t, u) : Perform one integration step.
        nstep(fun, u, t0, dt, nt) : Perform n times integration steps.

    References:
    .. [1] J.H Williamson, "Low-storage Runge-Kutta schemes", 1980.
    """
    nstage = 3
    butcherTable = np.array([[0, 0, 0, 0],
                             [1/3, 1/3, 0, 0],
                             [3/4, -3/16, 15/16, 0],
                             [0, 1/6, 3/10, 8/15]])
