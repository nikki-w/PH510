#!/usr/bin/env python3
"""MIT License

Copyright (c) [2025] [Nikki Walker]

Assignment 4: MC evaluation of Green's functions"""
import math
import numpy as np
from mpi4py import MPI
import montecarlo_class as MC

# Task 1:
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

        charge_dist: Funct that takes grid points and returns the charge
        """
        for i in range(self.n):
            for j in range(self.n):
                self.f[i, j] = charge_dist(i, j)

    def overrelaxation(self):
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

    def convergence_check(self, max_iter=10000, tolerance=1e-6):
        """Check for convergence"""
        for iteration in range(max_iter):
            init_phi = self.phi.copy()
            self.overrelaxation()
            delta = np.max(np.abs(self.phi - init_phi))
            # For convergance
            if delta < tolerance:
                return iteration #stop loop when delta falls below tolerance
        # For non-convergance
        return max_iter

# Task 2:
# Create new class to solve Green's eqtn
class GreensEqtn:
    """Create class to solve Green's function via use of random walkers"""
    def __init__(self, n, h, boundary_cond, charge_dist):
        """
        Initialise Green's solver

        n: grid size
        h: grid spacing
        boundary_cond: function that returns potential and boundary points
        charge_dist: function returning charge dist at grid points
        """
        self.n = n
        self.h = h
        self.boundary_cond = boundary_cond
        self.charge_dist = charge_dist
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()

    def random_walkers(self, start_i, start_j, n_walks):
        """
        Function that produces random walks from start points (i and j) to estimate Green's
        function

        grid_points_count: Counts how many walks pass through each grid point
        bound_end: Counts how many random walks will end at each boundary point
        """
        grid_points_count = np.zeros((self.n, self.n))
        bound_end = np.zeros((self.n, self.n))

        # Create random walks
        k = 0
        for k in range(n_walks):
            i, j = start_i, start_j
            while True:
                if i== 0 or i == (self.n-1)or j == 0 or j == (self.n-1):
                    bound_end[i, j] += 1
                    k += 1
                    break #stops random walk when it reaches the grid boundary

                # Randomly move to up, down, left or right
                position = np.random.randint(0, 4)
                if position == 0:
                    i += 1
                elif position == 1:
                    i -= 1
                elif position == 2:
                    j += 1
                else:
                    j -= 1

                # Sum visits through grid points
                if 0 < i < (self.n-1) and 0 < j < (self.n-1):
                    grid_points_count[i, j] += 1

        g_lap = bound_end / n_walks
        # From hints
        g_charge = (self.h**2 * grid_points_count) / n_walks

        return g_lap, g_charge

    def potential(self, start_i, start_j, n_walks=10000):
        """
        Calculates the potential using greens function
        """
        g_lap, g_charge = self.random_walkers(start_i, start_j, n_walks)

        # Laplace part of equation
        boundary_sum = 0
        for i in [0, (self.n-1)]:
            for j in range(self.n):
                boundary_sum += g_lap[i, j] * self.boundary_cond(i, j)

        for j in [0, (self.n-1)]:
            for i in range(1, (self.n-1)):
                boundary_sum += g_lap[i, j] * self.boundary_cond(i, j)

        # Poisson part of equation
        charge_sum = 0
        for i in range(1, (self.n-1)):
            for j in range(1, (self.n-1)):
                charge_sum += g_charge[i, j] * self.charge_dist(i, j)

        # Calculate total potential (sum laplace and poisson)
        phi = boundary_sum + charge_sum

        return phi
