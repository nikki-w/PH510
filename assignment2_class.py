#!/usr/bin/env python3
"""Python file that creates a class that can be 
imported to compute vector operations"""
import math

class Vector:
    """Object orientated vector creation"""
    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c

    def __str__(self):
        """Makes vectors printable"""
        return f"Vector: ({self.a:.2f}, {self.b:.2f}, {self.c:.2f})"

    def __add__(self, other):
        """Adds vector"""
        return Vector(self.a + other.a, self.b + other.b, self.c + other.c)

    def __sub__(self, other):
        """subtracts vector"""
        return Vector(self.a - other.a, self.b - other.b, self.c - other.c)

    def norm(self):
        """Calculates the magnitude of vector"""
        return math.sqrt(self.a**2 + self.b**2 + self.c**2)

    def dot(self, other):
        """Calculates the salar dot product"""
        return ((self.a * other.a) + (self.b * other.b) + (self.c * other.c))

    def cross(self, other):
        """Calculates the vector cross product of two vectors"""
        return Vector(((self.b * other.c) - (self.c * other.b)),
	((self.c * other.a) - (self.a * other.c)),
	((self.a * other.b) - (self.b * other.a)))
