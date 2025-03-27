#!/usr/bin/env python3
"""File conatining a Monte Carlo method defined as a class that is capable of
initialising independent parallel random number generators on workers from lead node
with an optional user defined global seed. The method is able to allow users to collect statistics"""
import math
import numpy as np
from mpi4py import MPI

class MC:
    def __init__(self, N, D, D_low=-1.0, D_high=1.0, seed=None):
        """
        Initialsie Monte Carlo (MC)
        N - Number of samples
        D - Number of dimensions for each point
        D_low - Lower limit for each dimension (set to -1.0 as default)
        D_high - Upper limit for each dimension (Set to 1.0 as default)
        """
        comm = MPI.COMM_WORLD
        self.comm = comm
        self.rank = comm.Get_rank()
        self.size = comm.Get_size()

        self.N = N
        self.D = D
        self.D_low = D_low
        self.D_high = D_high

        if seed is not None:
            self.seed = seed + self.rank
        else:
            self.seed = 42 + self.rank

        self.rng = np.random.default_rng(self.seed)

    def points(self):
        """Generates random points with D being the number of dimensions"""
        return self.rng.uniform(self.D_low, self.D_high, (self.N, self.D))

    


