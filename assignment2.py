#!/usr/bin/env python3
"""Python file that imports class with the use 
of performing operations with vectors"""
import assignment2_class as vect

print("\b")
print("Cartesian:")

v1 = vect.Vector(2, 10, 7)
print("First vector", v1)

v2 = vect.Vector(5, -2, 3)
print("Second vector", v2)

print("Sum of vectors:", v1 + v2)
print("Subtraction of vectors", v1 - v2)
print("Magnitude of first vector", v1.norm())
print("Dot product of two vectors", v1.dot(v2))
print("Cross product of two vectors", v1.cross(v2))

print("\b")
print("Spherical Polar:")
v1_sp = vect.SphericalPolar.cartesian_to_spherical(v1)
print("First vector:", v1_sp)

v2_sp = vect.SphericalPolar.cartesian_to_spherical(v2)
print("Second vector:", v2_sp)

print("Sum of vectors:", v1_sp + v2_sp)
print("Subtraction of vectors", v1_sp - v2_sp)
print("Magnitude of first vector", v1_sp.norm())
print("Dot product of two vectors", v1_sp.dot(v2_sp))
print("Cross product of two vectors", v1_sp.cross(v2_sp))
