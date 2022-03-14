import numpy as np

from integrate.lowStorageRungeKuttaMethods import LSRK3

class channel(object):
    def __init__(self, nu):
        self.nu = nu
        self.integrate = LSRK3()
