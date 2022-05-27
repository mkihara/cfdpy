import numpy as np


def wallNormalGrid(num:int, eta=1., delta=1.):
    """The non-uniformly spaced grid for the wall-normal direction `[1]`_.
    
    Args:
        num (int): Number of grids to generate.
        eta (float, optional): The parameter to controls how strongly the points are clusterd near the wall. Defaults to 1..
        delta (float, optional): The channel half-width. Defaults to 1..

    Returns:
        np.array: Non-uniformly spaced grid.

    .. _[1] [Lee and Moser, "Direct numerical simulation of turbulent channel flow up to Ret=5200", 2015.](https://arxiv.org/abs/1410.7809)
    """
    zeta = np.linspace(start=-1., stop=1., num=num)
    return delta * np.sin(eta*zeta*np.pi*0.5) / np.sin(eta*np.pi*0.5)


# def wallNormalGrid2(num:int, alpha=None):
#     """The non-uniformly spaced grid for the wall-normal direction [#2]_.

#     Args:
#         num (int): Number of grids to generate.
#         alpha (float): Adjustable patameter. 0.967 for Ret=180, 0.98 for Ret=395, None for Ret=640

#     Returns:
#         _type_: _description_
    
#     ToDo:
#         Something wrong with this method.

#     .. [#2]  Hiroyuki Abe, Hiroshi Kawamura, Yuichi Matsuo, "Direct Numerical Simulation of a Fully Developed Turbulent Channel Flow With Respect to the Reynolds Number Dependence", 2001.
#     """    
#     j = np.arange(num)
#     zeta = 2 * j / num - 1
#     if alpha==None:
#         alpha = 0.9885 - 0.5*zeta**2 + 0.405*zeta**3
#     return np.tanh(zeta / np.tanh(alpha)) / (2*alpha) + 0.5


if __name__=="__main__":
    Ret = 1000
    ny = 513
    eta = 0.99
    delta = 1.
    y = wallNormalGrid(ny, eta, delta)
    dy = np.diff(y)
    dyplus_wall = np.abs(dy).min() * Ret
    dyplus_center = np.abs(dy).max() * Ret
    print(dyplus_wall, dyplus_center)
    assert dyplus_wall < 0.45
    assert dyplus_center < 6.5
    # growthRate_wall = dy[1] / dy[0]
    # print(growthRate_wall)

    nx = 640
    Lx = 2*np.pi
    dx = Lx / nx
    dxplus = dx * Ret
    print(dxplus)
    assert dxplus < 10

    nz = 540
    Lz = np.pi
    dz = Lz / nz
    dzplus = dz * Ret
    print(dzplus)
    assert dzplus < 6
