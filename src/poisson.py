import numpy as np

from integrate.lowStorageRungeKuttaMethods import LSRK3

class poisson(object):
    def __init__(self):
        self.integrate = LSRK3()
