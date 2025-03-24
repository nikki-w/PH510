#!/usr/bin/env python3
"""File used for the creation of Monte Carlo class"""
import math
import numpy as np
from mpi4py import MPI

class MC:
    def __init__(self, N, D, D_low, D_high):
    """
    Initialsie Monte Carlo (MC)
        
    N - Number of samples
    D - Number of dimensions for each point
    D_low - Lower limit for each dimension
    D_high - Upper limit for each dimension
    """
    self.N = N
    self.D = D
    self.D_low = D_low
    self.D_high = D_high


