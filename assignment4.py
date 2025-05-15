#!/usr/bin/env python3
"""MIT License

Copyright (c) [2025] [Nikki Walker]

Assignment 4: MC evaluation of Green's functions"""
import math
import numpy as np
import montecarlo_class as MC

# Create new class to solve Poisson's equation
class PoissonEqtn:
    """ Create a class to solve Poissons equation for NxN grid via
        relaxation/over-relaxation
    """

    def __init__(self, grid_size, h, omega=None):
        """
        Initialises poissons eqtn solver

        grid_size: NxN grid
        h: Grid spacing (meters)
        omega: over-relaxation param
        """
        self.grid_size = grid_size
        self.h = h
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size

        # From hints, set an optimal omega when value isn't
        # defined
        if omega is None:
            self.omega = 2 / (1 + np.sin(math.pi/grid_size))
        else:
            self.omega = omega
