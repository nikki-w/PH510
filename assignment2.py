#!/usr/bin/env python3
"""Python file that imports class with the use 
of performing operations with vectors"""
import math
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
#v1_sp = vect.SphericalPolar(5, math.pi/3, math.pi/2)
#print("First vector:", v1_sp)

#v2_sp = vect.SphericalPolar(2, math.pi/5, math.pi/3)
#print("Second vector:", v2_sp)

#print("Sum of vectors:", v1_sp + v2_sp)



