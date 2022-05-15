import numpy as np

from cfdpy.integrate.lowStorageRungeKuttaMethods import LSRK3

class channel(object):
    """_summary_
    """    
    def __init__(self, nu):
        self.integrate = LSRK3
        raise NotImplementedError()
