import numpy as np
from scipy.sparse.linalg import LinearOperator, cg, gmres, lgmres

from derivative.finiteDifferenceMethods import centralFDM


class poisson2(object):
    def __init__(self, nx, ny, dx=1., dy=1.):
        self.nx = nx
        self.ny = ny
        self.dx = dx
        self.dy = dy
        self.p_shape = (self.ny, self.nx)
        self.A_shape = (self.nx*self.ny, self.nx*self.ny)
        self.derivative = centralFDM(order=6, highestDerivative=2)
        self.exitCode = 0
        self.A = LinearOperator(self.A_shape, matvec=self.lhs)

    def lhs(self, p_flatten):
        p = np.reshape(p_flatten, self.p_shape)
        p = ( self.derivative(p, h=self.dx, axis=-1, derivative=2) 
            + self.derivative(p, h=self.dy, axis=-2, derivative=2))
        return p.flatten()

    def rhs(self, u, v):
        return ( self.derivative(u, h=self.dx, axis=-1)
               + self.derivative(v, h=self.dy, axis=-2)).flatten()

    def solve(self, u, v, method="gmres"):
        if method=="gmres":
            p_flatten, self.exitCode = gmres(A=self.A, b=self.rhs(u, v))
        elif method=="cg":
            p_flatten, self.exitCode = cg(A=self.A, b=self.rhs(u, v))
        elif method=="lgmres":
            p_flatten, self.exitCode = lgmres(A=self.A, b=self.rhs(u, v))
        else:
            p_flatten = np.linalg.solve(self.A(np.identity(self.nx*self.ny)), self.rhs(u, v))
        return np.reshape(p_flatten, self.p_shape)
