"""Spectral methods
"""

import numpy as np
from math import ceil
from numpy import pi
from numpy.fft import fft, ifft, fftfreq, fftshift, ifftshift


class spectralMethod1(object):
    """Spectral Method.

    Args:
        n (int): The length of the array along the axis specified by axis.
        d (float): Sample spacing (inverse of the sampling rate). Defaults to 1.
        axis: Axis over which to compute the FFT. Defaults to 1.
        norm: Normalization mode. Default is "forward".
    
    Attributes:
        k (np.array): Wave Numbers.
        pad_width (int): Number of values padded to the edges of the specified axis.
    """    
    def __init__(self, n:int, d=1., axis=-1, norm='forward'):
        self.n = n
        self.d = d
        self.axis = axis
        self.norm = norm
        self.k = fftfreq(n, d) * 2*pi
        self.pad_width = ceil(self.n / 4)

    def fft(self, u):
        """Compute the discrete Fourier Transform.

        Args:
            u (np.array): State in physical space.

        Returns:
            np.array: State in wave space.
        """
        return fft(u, axis=self.axis, norm=self.norm)

    def ifft(self, u):
        """Compute the inverse discrete Fourier Transform.

        Args:
            u (np.array): State in physical space.

        Returns:
            np.array: State in wave space.
        """
        return ifft(u, axis=self.axis, norm=self.norm)

    def fftshift(self, u):
        """Shift the zero-frequency component to the center of the spectrum.

        Args:
            u (np.array): Input array.

        Returns:
            np.array: The shifted array.
        """
        return fftshift(u, self.axis)

    def ifftshift(self, u):
        """The inverse of fftshift.

        Args:
            u (np.array): Input array.

        Returns:
            np.array: The shifted array.
        """
        return ifftshift(u, self.axis)

    def diff(self, u, order:int=1):
        """Differentiate in the wave space.

        Args:
            u (np.array): States in the wave space.
            order (int): Order of the differential.

        Returns:
            np.array: Differentiated states in the wave space.
        """
        dim = np.ndim(u)
        k_shape = [1]*dim
        k_shape[self.axis] = len(self.k)
        k_shape = tuple(k_shape)
        k_expand = self.k.reshape(k_shape)
        return u * (1j * k_expand)**order

    def diff_phys(self, u, order:int=1):
        """Differentiate in the wave space.

        Args:
            u (np.array): States in the physical space.
            order (int): Order of the differential.

        Returns:
            np.array: Differentiated states in the physical space.
        """
        u = self.fft(u)
        u = self.diff(u, order)
        u = self.ifft(u)
        return u.real

    def multiply(self, u1, u2):
        """Multiply arguments element-wise.
        The aliasing error is eliminated by zero-padding (3/2 rule).

        Args:
            u1 (np.array): States in the wave space.
            u2 (np.array): States in the wave space.

        Returns:
            np.array: u1 * u2 in the wave space.
        """
        dim = np.ndim(u1)
        pad_width_ = np.zeros((dim,2), dtype=int)
        pad_width_[self.axis,:] = self.pad_width
        pad_width_ = tuple(map(tuple, pad_width_))

        ### padding
        u1 = self.fftshift(u1)
        u2 = self.fftshift(u2)
        p1 = np.pad(u1, pad_width_, mode='constant', constant_values=0.)
        p2 = np.pad(u2, pad_width_, mode='constant', constant_values=0.)
        p1 = self.ifftshift(p1)
        p2 = self.ifftshift(p2)

        ### transform states from wave space to 3/2 physical space
        p1 = self.ifft(p1)
        p2 = self.ifft(p2)

        ### multiply
        p1p2 = p1 * p2

        ### transform states from 3/2 physical space to 3/2 wave space
        p1p2 = self.fft(p1p2)

        ### unpadding
        p1p2 = self.fftshift(p1p2)
        mask = np.pad(np.ones_like(u1), pad_width_, mode='constant', constant_values=0.)
        p1p2 = np.compress(condition=mask, a=p1p2)
        p1p2 = self.ifftshift(p1p2)
        return p1p2

    def multiply_phys(self, u1, u2):
        """Multiply arguments element-wise.
        The aliasing error is eliminated by zero-padding (3/2 rule).

        Args:
            u1 (np.array): States in the physical space.
            u2 (np.array): States in the physical space.

        Returns:
            np.array: u1 * u2 in the physical space.
        """
        u1 = self.fft(u1)
        u2 = self.fft(u2)
        u1u2 = self.multiply(u1, u2)
        u1u2 = self.ifft(u1u2)
        return u1u2.real
