#!/usr/bin/env python3
"""Python file that creates x class that can be 
imported to compute vector operations"""
import math


# Create class for task 1
class Vector:
    """Object orientated vector creation"""
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        """Makes vectors printable"""
        return f"Vector: ({self.x:.2f}, {self.y:.2f}, {self.z:.2f})"

    def __add__(self, other):
        """Adds vector"""
        return Vector(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        """subtracts vector"""
        return Vector(self.x - other.x, self.y - other.y, self.z - other.z)

    def norm(self):
        """Calculates the magnitude of vector"""
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    def dot(self, other):
        """Calculates the salar dot product"""
        return ((self.x * other.x) + (self.y * other.y) + (self.z * other.z))

    def cross(self, other):
        """Calculates the vector cross product of two vectors"""
        return Vector(((self.y * other.z) - (self.z * other.y)),
	((self.z * other.x) - (self.x * other.z)),
	((self.x * other.y) - (self.y * other.x)))

    def cartesian_to_spherical(self):
        """Converts cartesian coordinates to spherical polar"""
        r = math.sqrt(self.x**2 + self.y**2 + self.z**2)
        theta = math.acos(self.z / r)
        phi = math.atan2(self.y, self.z)
        return r, theta, phi

    def spherical_to_cartesian(self):
        """Converts spherical polar coordinates to cartesian"""
        x = r * math.sin(theta) * math.cos(phi)
        y = r* math.sin(theta) * math.sin(phi)
        z = r * math.cos(theta)
        return x, y, z




