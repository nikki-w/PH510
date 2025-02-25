#!/usr/bin/env python3
"""Python file that creates x class that can be 
imported to compute vector operations"""
import math


# Create class for task 1
class Vector:
    """Creates vectors in cartesian coordinates and performs operations"""
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
        """Subtracts vector"""
        return Vector(self.x - other.x, self.y - other.y, self.z - other.z)

    def norm(self):
        """Calculates the magnitude of vector"""
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    def dot(self, other):
        """Calculates the scalar dot product"""
        return ((self.x * other.x) + (self.y * other.y) + (self.z * other.z))

    def cross(self, other):
        """Calculates the vector cross product of two vectors"""
        return Vector(((self.y * other.z) - (self.z * other.y)),
	((self.z * other.x) - (self.x * other.z)),
	((self.x * other.y) - (self.y * other.x)))

    def cartesian_to_spherical(self):
        """Converts cartesian coordinates to spherical polar coordinates"""
        r = self.norm()
        theta = math.acos(self.z / r) if r != 0 else 0
        phi = math.atan(self.y / self.x) if self.x != 0 else 0
        return SphericalPolar(r, theta, phi)

class SphericalPolar(Vector):
    """Inherits Vector class and converts into Spherical polar coordinates"""
    def __init__(self, r, theta, phi):
        """Initialise"""
        self.r = r
        self.theta = theta
        self.phi = phi
        # Call __init__ from Vector
        super().__init__(self.r, self.theta, self.phi)

    def __str__(self):
        """Makes vector printable"""
        return f"Spherical Polar: ({self.r:.2f}, {self.theta:.2f}, {self.phi:.2f})"

    def spherical_to_cartesian(self):
        """Converts spherical polar to cartesian coordinates"""
        x = self.r * math.sin(self.theta) * math.cos(self.phi)
        y = self.r * math.sin(self.theta) * math.sin(self.phi)
        z = self.r * math.cos(self.theta)
        return Vector(x, y, z)

    def __add__(self, other):
        """Adds two spherical polar vectors"""
        add_sp = self.spherical_to_cartesian() + other.spherical_to_cartesian()
        return add_sp.cartesian_to_spherical()

    def __sub__(self, other):
        """Adds two spherical polar vectors"""
        sub_sp = self.spherical_to_cartesian() - other.spherical_to_cartesian()
        return sub_sp.cartesian_to_spherical()

    def norm(self):
        """Calculates the magnitude of vector"""
        return self.r

    def dot(self, other):
        """Calculates the salar dot product"""
        dot_sp = self.spherical_to_cartesian().dot(other.spherical_to_cartesian())
        return dot_sp

    def cross(self, other):
        """Calculates the salar dot product"""
        cross_sp = self.spherical_to_cartesian().cross(other.spherical_to_cartesian())
        return cross_sp.cartesian_to_spherical()
