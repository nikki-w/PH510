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

        # Create arrays to store potential, phi, and charge, f, for later tasks
        self.phi = np.zeros((n, n))
        self.f = np.zeros((n, n))

        # Split up tasks for parallel implementation
        self.localrow = n // self.size
        self.extrarows = n % self.size
        self.startrow = self. startrow + self.localrow + (1 if self.rank < self.extrarows else 0)

    def bc_potential(self, boundary_cond):
        """
        Define boundary conditions to place the potential

        boundary_cond: funct that takes grid points and returns the potential
        """
        # Define boundary values
        # Left and right boundary:
        for j in [0, (self.n - 1)]:
            for i in range(self.n):
                self.phi[i, j] = boundary_cond(i, j)
        # Top and bottom boundary:
        for i in [0, (self.n - 1)]:
            for j in range(self.n):
                self.phi[i, j] = boundary_cond(i, j)

    def charge_dist(self, charge_dist):
        """
        Define the distribution of charge across the grid

        charge_dist: Funct that takes grid points and returnes the charge
        """
        for i in range(self.n):
            for j in range(self.n):
                self.f[i, j] = charge_dist(i, j)

    def overelaxation(self):
        """
        Method for relaxation/overrelaxation
        """
        phi_new = self.phi.copy() # copy phi values into new variable
        # Method will update the potential via overrelaxation for points inside grid
        for i in range(1, (self.n - 1)):
            for j in range(1, (self.n - 1)):
                # Define new phi from equation given in hand-out
                terms = (self.phi[i+1, j] + self.phi[i-1, j] + self.phi[i, j+1] + self.phi[i, j-1])
                charge = terms + self.h**2 * self.f[i, j]
                phi_new = ((1 - self.omega) * self.phi[i, j] + (self.omega/4) * charge)
        self.phi = phi_new


