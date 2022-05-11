import numpy as np

from integrate.lowStorageRungeKuttaMethods import LSRK3

class channel(object):
    """_summary_
    """    
    def __init__(self, nu):
        self.integrate = LSRK3
        raise NotImplementedError()
