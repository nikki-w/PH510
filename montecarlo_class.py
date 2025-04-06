#!/usr/bin/env python3
"""MIT License

Copyright (c) [2025] [Nikki Walker]

Importable module with a parallel Monte Carlo class allowing users to 
create simulations with parallel RNG, collect statistics, estimate volumes 
of n-dimensional spheres and integrate gaussian functions in n-dimensions"""

import math
import numpy as np
from mpi4py import MPI
from scipy.special import gamma

class MC:
    """Parallel MC class

        comm: MPI communicator
        rank: Process rank
        size: Number of processes
        n: Number of samples
        d: Number of dimensions
        d_low: Lower limit (integration)
        d_high: Upper limit (integration)
        seed: Random seed
        rng: Random number generator
    """

    def __init__(self, n, d, d_low=-1.0, d_high=1.0, seed=None):
        """Initialsie parallel Monte Carlo (MC).

        n: Number of samples across all processes
        d: Number of dimensions for each point
        d_low: Lower limit for each dimension (set to -1.0 as default)
        d_high: Upper limit for each dimension (set to 1.0 as default)
        seed: Optional global seed for RNG (set to None as default)
        """
        comm = MPI.COMM_WORLD
        self.comm = comm
        self.rank = comm.Get_rank()
        self.size = comm.Get_size()

        # Utilise MPI to distribute samples evenly across processes
        self.n = n
        self.n_local = n // self.size
        if self.rank < n % self.size:
            self.n_local += 1

        self.d = d
        self.d_low = d_low
        self.d_high = d_high

        #Initialise rng using seed
        if seed is not None:
            self.seed = (seed + self.rank) 
        else:
            self.seed = None

        self.rng = np.random.default_rng(self.seed)

    def points(self):
        """Generates random points in d-dimensions

        Returns:
            np.darray: Array of shape (n_local, d) with random points
        """
        return self.rng.uniform(self.d_low, self.d_high, (self.n_local, self.d))

    def volume(self, r):
        """Estimates the volume of a unit ball in d-dimensions with radius, r

        Params:
            r: Radius of sphere (set to 1.0 as default)

        Returns:
            tuple: Volume estimate
        """
        dist_squared = np.sum(self.points()**2, axis=1)
        points_inside = np.sum(dist_squared <= r**2)

        # Calculate the estimated volume
        volume_est = (points_inside/self.n_local) * ((self.d_high - self.d_low)**self.d)
        return volume_est

    def volume_exact(self, r=1):
        """Calculates the exact volume of a unit ball in d dimensions

        Params:
            r: Raduis of sphere (set to 1.0 as default)

        Returns:
            float: Exact Volume
        """
        return (math.pi**(self.d/2) / gamma((self.d/2) + 1)) * r**self.d

    def gauss_integ(self, sigma, x0):
        """ doc for integral"""
        sigma = np.asarray(sigma)
        x0 = np.asarray(x0)

        # Have error raised if sigma and x0 are not equal to the number of dimensions
        if len(sigma) != self.d or len(x0) != self.d:
            raise ValueError("Dimensions of both sigma and x0 need to equal the number of dimensions, d")

        t = self.points()
        x = t / (1-t**2)
        a = np.prod((1+t**2) / (1-t**2)**2, axis=1)
        exponent = -0.5 * np.sum((x-x0)**2 / sigma**2, axis=1)
        gauss = np.exp(exponent) / ((2*np.pi)**(self.d/2) * np.prod(sigma))

        integ = gauss * a
        return self.statistics(integ, is_integral=True)

    def statistics(self, local_value, is_integral=False):
        """Collects statistics such as mean and variance across MPI processes

        Params:
            local_values: local values to collect statistics for

        Returns:
            tuple: (mean, error) if rank is 0, else (None, None)
        """
        if is_integral:
            local_sum = np.sum(local_value)
            local_sum_squ = np.sum(local_value**2)
            normalisation = (self.d_high - self.d_low)**self.d
        else:
            local_sum = np.sum(local_value)
            local_sum_squ = np.sum(local_value**2)
            normalisation = 1.0

        # Reduce across all MPI processws
        global_sum = self.comm.reduce(local_sum, op=MPI.SUM, root=0)
        global_sum_squ = self.comm.reduce(local_sum_squ, op=MPI.SUM, root=0)

        if self.rank == 0:
            mean = (global_sum / self.n) * normalisation
            variance = (global_sum_squ / self.n) * normalisation**2 - mean**2
            error = math.sqrt(variance / self.n)
            return mean, error, variance
        return None, None, None
