"""Finite Difference Methods
"""

import numpy as np


def FDMWeights(M, x0, alpha):
    """Calculate the weights in finite difference formulas
    for any order of derivative and to any order of accuracy
    on onedimensional grids with arbitrary spacing.

    Args:
        M (int): Order of derivative
        x0 (float): Approximations at this point
        alpha (np.array): x-cordinates. length must be N
    
    Attributes:
        N (int): Order of accuracy, which is equivalent to len(alpha)-1.
    
    Returns:
        np.array: Weights

    References: 
        Bengt Fornberg, "Generation of Finite Difference Formulas on Arbitrarily Spaced Grids", 1988.
    """
    N = len(alpha) - 1
    delta = np.zeros([M+1,N+1,N+1])
    delta[0,0,0] = 1.
    c1 = 1.
    for n in range(1, N+1):
        c2 = 1.
        for nu in range(n):
            c3 = alpha[n] - alpha[nu]
            c2 *= c3
            for m in range(min(n, M)+1):
                delta[m,n,nu] = ((alpha[n]-x0)*delta[m,n-1,nu] - m*delta[m-1,n-1,nu]) / c3
        for m in range(min(n, M)+1):
            delta[m,n,n] = c1/c2 * (m*delta[m-1,n-1,n-1] - (alpha[n-1]-x0)*delta[m,n-1,n-1])
        c1 = c2
    return delta


class centralFDM(object):
    """Central Finite Difference Method

    Args:
        order (int, optional): The order of the accuracy. Defaults to 2.
        highestDerivative (int, optional): The order of the highest derivative. Defaults to 1.
    """    
    def __init__(self, order:int=2, highestDerivative=1):  
        assert (order % 2) == 0, "order must be even number."
        assert order > 0, "order must be greater than 0."
        assert highestDerivative > 0, "highestDerivative must be greater than 0."
        self.order = order
        self.highestDerivative = highestDerivative
        self.nGridPoints = ((self.highestDerivative + 1) // 2) * 2 - 1 + self.order
        self.set_alpha()
        self.weight = FDMWeights(M=self.highestDerivative, x0=0, alpha=self.alpha)[:,self.order]

    def __call__(self, f, axis=-1, derivative=1, h=1.):
        """Calculate the derivative.

        Args:
            f (np.array): An array containing samples.
            axis (int, optional): The derivative is calculated only along the given axis. Defaults to -1.
            derivative (int, optional): The order of the derivative. Defaults to 1.
            h (float, optional): The space of the uniform grid. Defaults to 1..

        Returns:
            np.array: The derivative.
        """        
        df = np.zeros_like(f)
        weight_ = self.weight[derivative]
        alpha_ = self.alpha[weight_!=0]
        weight_ = weight_[weight_!=0]
        for i, alpha_i in enumerate(alpha_):
            df += np.roll(f, shift=-int(alpha_i), axis=axis) * weight_[i]
        return df / h**derivative

    def set_alpha(self):
        alpha_ = np.arange(self.nGridPoints, dtype=float)
        alpha_ = self.__infiniteSeries(alpha_)
        self.alpha = np.cumsum(alpha_)

    def __infiniteSeries(self, n):
        return n * (-1)**(n-1)


class upwindFDM(object):
    """Upwind Finite Difference Method

    Args:
        order (int, optional): The order of the accuracy. Defaults to 1.
        highestDerivative (int, optional): The order of the highest derivative. Defaults to 1.
    """  
    def __init__(self, order:int=1, highestDerivative:int=1):
        assert order > 0, "order must be greater than 0."
        assert highestDerivative > 0, "highestDerivative must be greater than 0."
        self.order = order
        self.highestDerivative = highestDerivative
        self.nGridPoints = self.order+self.highestDerivative
        self.start = - (self.nGridPoints) // 2
        self.alpha = np.arange(start=self.start, stop=self.start+self.nGridPoints)
        self.weight = FDMWeights(M=self.highestDerivative, x0=0., alpha=self.alpha)[:,self.order]
        self.weight2 = FDMWeights(M=self.highestDerivative, x0=0., alpha=-self.alpha)[:,self.order]

    def __call__(self, f, axis=-1, derivative=1, h=1., c=None):
        """Calculate the derivative.

        Args:
            f (np.array): An array containing samples.
            axis (int, optional): The derivative is calculated only along the given axis. Defaults to -1.
            derivative (int, optional): The order of the derivative. Defaults to 1.
            h (float, optional): The space of the uniform grid. Defaults to 1..
            c (float or np.array, optional): The advection speed. Defaults to None.

        Returns:
            np.array: The derivative.
        """
        df = np.zeros_like(f)
        df2 = np.zeros_like(f)
        for i, alpha_i in enumerate(self.alpha):
            df += np.roll(f, shift=-int(alpha_i), axis=axis) * self.weight[derivative,i]
            df2 += np.roll(f, shift=int(alpha_i), axis=axis) * self.weight2[derivative,i]
        if c == None:
            c = f
        df = np.where(c>=0, df, df2)
        return df / h**derivative
