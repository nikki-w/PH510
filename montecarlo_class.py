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
        #points = self.points()
        dist_squared = np.sum(self.points()**2, axis=1)
        points_inside = np.sum(dist_squared <= r**2)

        if self.rank == 0:
            print(f"Points inside: {points_inside}/{self.n_local}")
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

    def statistics(self, local_value):
        """Collects statistics such as mean and variance across MPI processes

        Params:
            local_values: local values to collect statistics for

        Returns:
            tuple: (mean, error) if rank is 0, else (None, None)
        """
        local_sum = np.sum(local_value)
        local_sum_squ = np.sum(local_value**2)

        # Reduce across all MPI processws
        global_sum = self.comm.reduce(local_sum, op=MPI.SUM, root=0)
        global_sum_squ = self.comm.reduce(local_sum_squ, op=MPI.SUM, root=0)

        if self.rank == 0:
            mean = global_sum / self.n
            variance = (global_sum_squ / self.n) - mean**2
            error = math.sqrt(variance / self.n)
            return mean, error
        return None, None
