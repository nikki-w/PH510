#!/usr/bin/env python3
"""Python3 script using mpi4py to do mid point quadrature 
calculations to estimate a value of pi"""
import time
from mpi4py import MPI

# Start time
time0 = time.time()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
nprocs = comm.Get_size()

N = 3000000000
DELTA = 1.0 / N

# Define integral
def integral(x):
    """Defininte integral used to estimate a value of pi 
    using mid point quadrature for values of x"""
    return 4.0 / (1.0 + x*x)

if rank == 0:
    # Want to split the caluclation of mpr evenly into sub-tasks, this is done below
    # by dividing N into nprocs (number of steps into number of processes)
    quo, rem = divmod(N, nprocs)
    # Distribute sub-tasks tasks evenly
    subs = []
    for i in range(nprocs):
        if i < rem:
            subs.append(quo + 1)
        else:
            subs.append(quo)

    # Define start and end points for integral
    start_integ = [sum(subs[:i]) for i in range(nprocs)]
    end_integ = [sum(subs[:i + 1]) for i in range(nprocs)]

    INTEG = [(start_integ[i], end_integ[i]) for i in range(nprocs)]
else:
    INTEG = None

# Scatter chunks of the calculation to different processes
INTEG = comm.scatter(INTEG, root=0)

# Compute different sections of pi across processes
PI = 0.0
for j in range(INTEG[0], INTEG[1]):
    mpr = (j + 0.5) * DELTA
    y = integral(mpr)
    PI += y

PI = PI * DELTA

# Gather and sum the data using comm.Reduce to form one complete pi estimate
PI = comm.reduce(PI, op=MPI.SUM, root=0)

if rank == 0:
    print(f'pi computed in {(time.time() - time0):.3f} sec.')
    print(f'Pi is {PI:.15}.')
