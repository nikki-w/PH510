import assignment2_class as vect
import math

v1 = vect.Vector(2, 10, 7)
print("First vector", v1)

v2 = vect.Vector(5, -2, 3)
print("Second vector", v2)

print("Sum of vectors:", v1 + v2)
print("Subtraction of vectors", v1 - v2)
print("Magnitude of first vector", v1.norm())
print("Dot product of two vectors", v1.dot(v2))
print("Cross product of two vectors", v1.cross(v2))
