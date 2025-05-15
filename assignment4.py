#!/usr/bin/env python3
"""MIT License

Copyright (c) [2025] [Nikki Walker]

Assignment 4: MC evaluation of Green's functions"""
import math
import numpy as np
from mpi4py import MPI
import montecarlo_class as MC

# Create new class to solve Poisson's equation
class PoissonEqtn:
    """ Create a class to solve Poissons equation for NxN grid via
        relaxation/over-relaxation
    """

    def __init__(self, n, h, omega=None):
        """
        Initialises poissons eqtn solver

        n: NxN grid
        h: Grid spacing (meters)
        omega: over-relaxation param
        """
        self.n = n
        self.h = h
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size

        # From hints, set an optimal omega when value isn't
        # defined
        if omega is None:
            self.omega = 2 / (1 + np.sin(math.pi/n))
        else:
            self.omega = omega

        # Create arrays to store potential, phi, and charge, fm for later tasks
        self.phi = np.zeros((n, n))
        self.f = np.zeros((n, n))

        # Split up tasks for parallel implementation
        self.localrow = n // self.size
        self.extrarows = n % self.size
        self.startrow = self. startrow + self.localrow + (1 if self.rank < self.extrarows else 0)
