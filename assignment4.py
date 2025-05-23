#!/usr/bin/env python3
"""MIT License

Copyright (c) [2025] [Nikki Walker]

Assignment 4: MC evaluation of Green's functions"""
import math
import numpy as np
import matplotlib.pyplot as plt
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
        self.size = self.comm.Get_size()

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
        self.extrarow = n % self.size
        self.startrow = self.rank * self.localrow
        if self.rank < self.extrarow:
            self.startrow += self.rank
            self.localrow += 1
        else:
            self.startrow += self.extrarow

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
        phi_new = np.copy(self.phi) # copy phi values into new variable

        # Method will update the potential via overrelaxation for points inside grid
        for i in range(1, self.n-1):
            for j in range(1, self.n-1):
                # Calculate avg
                phi_sum = self.phi[i+1, j] + self.phi[i-1, j] + self.phi[i, j+1] + self.phi[i, j-1]
                charge_term = min(self.h**2 * self.f[i, j], 1e10)
                update = (phi_sum + charge_term)*0.25
                delta = update - self.phi[i, j]
                phi_new[i, j] = self.phi[i, j] + min(self.omega * delta, 1e4)

                # Force finite values
                if not np.isfinite(phi_new[i,j]):
                    phi_new[i,j] = phi_sum

        self.phi = phi_new

    def convergence_check(self, max_iter=10000, tolerance=1e-6):
        """Check for convergence"""
        for iteration in range(max_iter):
            init_phi = self.phi.copy()
            self.overrelaxation()

            delta = np.nanmax(np.abs(self.phi - init_phi))
            # For convergance
            if np.isnan(delta) or delta < tolerance:
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
            steps = 0
            max_steps = 100 *self.n

            while steps < max_steps:
                if i== 0 or i == (self.n-1) or j == 0 or j == (self.n-1):
                    bound_end[i, j] += 1
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
                steps += 1

        # Calculate greens funct
        g_lap = bound_end / n_walks
        g_lap_err = np.sqrt(g_lap * (1-g_lap) / n_walks)

        # From hints
        g_charge = (self.h**2 * grid_points_count) / n_walks
        g_charge_err = np.sqrt(g_charge * (1 - g_charge) / n_walks)

        return g_lap, g_charge, g_lap_err, g_charge_err

    def potential(self, start_i, start_j, n_walks=10000):
        """
        Calculates the potential using greens function
        """
        g_lap, g_charge, g_lap_err, g_charge_err = self.random_walkers(start_i, start_j, n_walks)

        # Laplace part of equation
        boundary_sum = 0
        boundary_err = 0
        for i in [0, (self.n-1)]:
            for j in range(self.n):
                boundary_sum += g_lap[i, j] * self.boundary_cond(i, j)
                boundary_err += (g_lap_err[i, j] * self.boundary_cond(i, j))**2

        for j in [0, (self.n-1)]:
            for i in range(1, (self.n-1)):
                boundary_sum += g_lap[i, j] * self.boundary_cond(i, j)
                boundary_err += (g_lap_err[i, j] * self.boundary_cond(i, j))**2

        # Poisson part of equation
        charge_sum = 0
        charge_err = 0
        for i in range(1, (self.n-1)):
            for j in range(1, (self.n-1)):
                charge_sum += g_charge[i, j] * self.charge_dist(i, j)
                charge_err += (g_charge_err[i, j] * self.boundary_cond(i, j))**2

        # Calculate total potential (sum laplace and poisson)
        phi = boundary_sum + charge_sum
        err = np.sqrt(boundary_err + charge_err)

        return phi, err

    def potential_calc(self, g_lap, g_charge, g_lap_err=None, g_charge_err=None):
        """
        Calculates the potential using greens function (task 4)
        """
        if g_lap_err is None:
            g_lap_err = np.zeros_like(g_lap)
        if g_charge_err is None:
            g_charge_err = np.zeros_like(g_charge)

        # Laplace part of equation
        boundary_sum = 0
        boundary_err = 0
        for i in [0, (self.n-1)]:
            for j in range(self.n):
                boundary_sum += g_lap[i, j] * self.boundary_cond(i, j)
                boundary_err += (g_lap_err[i, j] * self.boundary_cond(i, j))**2

        for j in [0, (self.n-1)]:
            for i in range(1, (self.n-1)):
                boundary_sum += g_lap[i, j] * self.boundary_cond(i, j)
                boundary_err += (g_lap_err[i, j] * self.boundary_cond(i, j))**2

        # Poisson part of equation
        charge_sum = 0
        charge_err = 0
        for i in range(1, (self.n-1)):
            for j in range(1, (self.n-1)):
                charge_sum += g_charge[i, j] * self.charge_dist(i, j)
                charge_err += (g_charge_err[i, j] * self.boundary_cond(i, j))**2

        # Calculate total potential (sum laplace and poisson)
        phi = boundary_sum + charge_sum
        err = np.sqrt(boundary_err + charge_err)

        return phi, err

    def parallel_potential(self, start_i, start_j, tot_walks=100000):
        """Parallel computation of potential using MC class"""
        # Divide up walks across processes
        mc = MC.MC(tot_walks, d=1)
        n_walks_lo = mc.n_local

        # Generate walks for g_lap and g_charge
        g_lap, g_charge, g_lap_err, g_charge_err = self.random_walkers(start_i, start_j, n_walks_lo)

        # Reduce results across all processes
        global_g_lap = np.zeros((self.n, self.n))
        global_g_charge = np.zeros((self.n, self.n))
        global_g_lap_err = np.zeros((self.n, self.n))
        global_g_charge_err = np.zeros((self.n, self.n))
        self.comm.Reduce(g_lap, global_g_lap, op=MPI.SUM, root=0)
        self.comm.Reduce(g_charge, global_g_charge, op=MPI.SUM, root=0)
        self.comm.Reduce(g_lap_err**2, global_g_lap_err, op=MPI.SUM, root=0)
        self.comm.Reduce(g_charge_err**2, global_g_charge_err, op=MPI.SUM, root=0)

        # Normalise and generate results
        if mc.rank == 0:
            global_g_lap /= tot_walks
            global_g_charge /= tot_walks
            global_g_lap_err = np.sqrt(global_g_lap_err)
            global_g_charge_err = np.sqrt(global_g_charge_err)

            phi, err = self.potential_calc(global_g_lap, global_g_charge, global_g_lap_err, global_g_charge_err)
            return phi, err
        return None, None


# Task 3:
# Define parameters used in all parts of task
side_length = 0.1 #m
h = 0.001
n = int(side_length / h)

# Assume that the boundary condition is +1V for all of task 3
def bc_task3(i, j):
    if i == 0 or i == (n-1) or j == 0 or j == (n-1):
        return 1.0
    else:
        return 0

greens = GreensEqtn(n, h, bc_task3, lambda i, j: 0) # lambda funct to temp set i and j to 0

# Define list of points needed for a) - d)
points = [(n//2, n//2), (n//4, n//4), (1, n//4), (1, 1)]

if greens.rank == 0:
    print("Task 3: Green's Function Evaluation")

# Step through each point to determine answer
for k, (i, j) in enumerate(points):
    g_lap, g_charge, g_lap_err, g_charge_err = greens.random_walkers(i, j, 10000)

    if greens.rank == 0:
        # Results
        print(f"\nPoint {k+1} ({(i*h*100):.1f}cm, {(j*h*100):.1f}cm):")
        print(f"Max Laplace G: {np.max(g_lap):.4f} +/- {np.max(g_lap_err):.4f}")
        print(f"Max Charge G: {np.max(g_charge):.4f} +/- {np.max(g_charge_err):.4f}")

if greens.rank == 0:
    # Plot of results
    plt.figure()
    plt.imshow(g_lap, extent=[0, 10, 0, 10], origin='lower')
    plt.colorbar(label="Laplace Green's function")
    plt.title(f"Laplace Greens at ({(i*h*100):.1f}cm,{(j*h*100):.1f}cm)")
    plt.savefig(f"Task3_laplace_plot{k+1}.png")
    plt.close()

    plt.figure()
    plt.imshow(g_charge, extent=[0, 10, 0, 10], origin='lower')
    plt.colorbar(label="Charge Green's function")
    plt.title(f"Charge Greens at ({(i*h*100):.1f}cm,{(j*h*100):.1f}cm)")
    plt.savefig(f"Task3_charge_plot{k+1}.png")
    plt.close()

# Task 4:
# Parameters are the same as in Task 3 so no need to redefine
# Define new boundary conditions as per task requirements
def bc_task4_a(i, j):
    return 1.0

def bc_task4_b(i, j):
    if j == 0 or j == (n-1):
        return 1.0
    else:
        return -1.0

def bc_task4_c(i, j):
    if i == 0:
        return 2.0
    elif i == (n-1):
        return -4.0
    elif j == (n-1):
        return 2.0
    else:
        return 0

# Define charge distributions
def charge_task4_a(i, j):
    return 10/(n*n)

def charge_task4_b(i, j):
    return (n-1-j)/ (n-1)

def charge_task4_c(i, j):
    r = np.sqrt((i-n//2)**2 + ((j-n//2)**2)*h)
    return np.exp(-2000*r)

point_names = ["Centre", "Quarter", "Edge", "Corner"]

# No need to create new points as we use the ones from task 3
# Initialise poissons eqtn
poisson = PoissonEqtn(n, h, omega=1.7)
# Greens is already initialised above but need to redo so we can set the bc
greens = GreensEqtn(n, h, None, None)

if greens.rank == 0:
    print("\nTask 4: Potential Calculations")

# a)
if greens.rank == 0:
    print("\na): All edges +1V, uniform 10C charge")
greens.boundary_cond = bc_task4_a
greens.charge_dist = charge_task4_a
poisson.bc_potential(bc_task4_a)
poisson.charge_dist(charge_task4_a)
poisson.convergence_check()

for (i, j), name in zip(points, point_names):
    phi_mc, err_mc = greens.parallel_potential(i, j, 100000)
    if greens.rank == 0:
        print(f"{name} ({(i*h*100):.1f}cm,{(j*h*100):.1f}cm): "
                  f"MC = {phi_mc:.4f} +/- {err_mc:.4f} V, "
                  f"Relaxation = {poisson.phi[i,j]:.4f} V")


# b)
if greens.rank == 0:
    print("\nb): Top/bottom +1V, sides -1V, vertical gradiant")
greens.boundary_cond = bc_task4_b
greens.charge_dist = charge_task4_b
poisson.bc_potential(bc_task4_b)
poisson.charge_dist(charge_task4_b)
poisson.convergence_check()

for (i, j), name in zip(points, point_names):
    phi_mc, err_mc = greens.parallel_potential(i, j, 100000)
    if greens.rank == 0:
        print(f"{name} ({(i*h*100):.1f}cm,{(j*h*100):.1f}cm): "
                  f"MC = {phi_mc:.4f} +/- {err_mc:.4f} V, "
                  f"Relaxation = {poisson.phi[i,j]:.4f} V")

# c)
if greens.rank == 0:
    print("\nc): Mixed boundaries with exponential decay charge")
greens.boundary_cond = bc_task4_c
greens.charge_dist = charge_task4_c
poisson.bc_potential(bc_task4_c)
poisson.charge_dist(charge_task4_c)
poisson.convergence_check()

for (i, j), name in zip(points, point_names):
    phi_mc, err_mc = greens.parallel_potential(i, j, 100000)
    if greens.rank == 0:
        print(f"{name} ({(i*h*100):.1f}cm,{(j*h*100):.1f}cm): "
                  f"MC = {phi_mc:.4f} +/- {err_mc:.4f} V, "
                  f"Relaxation = {poisson.phi[i,j]:.4f} V")

# Task 5
if greens.rank == 0:
    print("\nTask 5: Validation")

    poisson_values = PoissonEqtn(n, h, omega=1.7)
    greens_values = GreensEqtn(n, h, None, None)

    #4a
    print("\n Task 4 a) validation:")
    poisson_values.bc_potential(bc_task4_a)
    poisson_values.charge_dist(charge_task4_a)
    poisson_values.convergence_check(max_iter=5000)

    greens_values.boundary_cond = bc_task4_a
    greens_values.charge_dist = charge_task4_a

    diff = []
    for (i, j), name in zip(points, point_names):
        phi_mc, err_mc = greens_values.parallel_potential(i, j, 100000)
        phi_relaxation = poisson_values.phi[i,j]
        diffs = abs(phi.mc - phi_relaxation)
        diff.append(diffs)
        print(f"{name}: (MC - Relax) = {diff:.2e} V")
