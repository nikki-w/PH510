#!/usr/bin/env python3
"""Python file that imports class with the use 
of performing operations with vectors"""
import math
import assignment2_class as vect

print("\b")
print("Task 1: Cartesian:")

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
print("Task 2: Spherical Polar:")
v1_sp = vect.SphericalPolar.cartesian_to_spherical(v1)
print("First vector:", v1_sp)

v2_sp = vect.SphericalPolar.cartesian_to_spherical(v2)
print("Second vector:", v2_sp)

print("Sum of vectors:", v1_sp + v2_sp)
print("Subtraction of vectors", v1_sp - v2_sp)
print("Magnitude of first vector", v1_sp.norm())
print("Dot product of two vectors", v1_sp.dot(v2_sp))
print("Cross product of two vectors", v1_sp.cross(v2_sp))

print("\b")
print("Task 3:")
print("a)")
triangle1_a = vect.Vector(0, 0, 0)
triangle1_b = vect.Vector(1, 0, 0)
triangle1_c = vect.Vector(0, 1, 0)
triangle1_ab = triangle1_b - triangle1_a
triangle1_ac = triangle1_c - triangle1_a
triangle1_ba = triangle1_a - triangle1_b
triangle1_bc = triangle1_c - triangle1_b
triangle1_ca = triangle1_a - triangle1_c
triangle1_cb = triangle1_b - triangle1_c
triangle1_area = 0.5 * (triangle1_ab.cross(triangle1_ac)).norm()
print("Area of Triangle 1 =", triangle1_area)

triangle2_a = vect.Vector(-1, -1, -1)
triangle2_b = vect.Vector(0, -1, -1)
triangle2_c = vect.Vector(-1, 0, -1)
triangle2_ab = triangle2_b - triangle2_a
triangle2_ac = triangle2_c - triangle2_a
triangle2_ba = triangle2_a - triangle2_b
triangle2_bc = triangle2_c - triangle2_b
triangle2_ca = triangle2_a - triangle2_c
triangle2_cb = triangle2_b - triangle2_c
triangle2_area = 0.5 * (triangle2_ab.cross(triangle2_ac)).norm()
print("Area of Triangle 2 =", triangle2_area)

triangle3_a = vect.Vector(1, 0, 0)
triangle3_b = vect.Vector(0, 0, 1)
triangle3_c = vect.Vector(0, 0, 0)
triangle3_ab = triangle3_b - triangle3_a
triangle3_ac = triangle3_c - triangle3_a
triangle3_ba = triangle3_a - triangle3_b
triangle3_bc = triangle3_c - triangle3_b
triangle3_ca = triangle3_a - triangle3_c
triangle3_cb = triangle3_b - triangle3_c
triangle3_area = 0.5 * (triangle3_ab.cross(triangle3_ac)).norm()
print("Area of Triangle 3 =", triangle3_area)

triangle4_a = vect.Vector(0, 0, 0)
triangle4_b = vect.Vector(1, -1, 0)
triangle4_c = vect.Vector(0, 0, 1)
triangle4_ab = triangle4_b - triangle4_a
triangle4_ac = triangle4_c - triangle4_a
triangle4_ba = triangle4_a - triangle4_b
triangle4_bc = triangle4_c - triangle4_b
triangle4_ca = triangle4_a - triangle4_c
triangle4_cb = triangle4_b - triangle4_c
triangle4_area = 0.5 * (triangle4_ab.cross(triangle4_ac)).norm()
print("Area of Triangle 4 =", triangle4_area)

print("b)")
# Triangle 1
# Angle AB and AC
triangle1_angle_a = math.acos((triangle1_ab.dot(triangle1_ac)) /
(triangle1_ab.norm() * triangle1_ac.norm()))
print("Triangle 1: Angle between AB and AC =", math.degrees(triangle1_angle_a))
# Angle BA and BC
triangle1_angle_b = math.acos((triangle1_ba.dot(triangle1_bc)) /
(triangle1_ba.norm() * triangle1_bc.norm()))
print("Triangle 1: Angle between BA and BC =", math.degrees(triangle1_angle_b))
# Angle CA and CB
triangle1_angle_c = math.acos((triangle1_ca.dot(triangle1_cb)) /
(triangle1_ca.norm() * triangle1_cb.norm()))
print("Triangle 1: Angle between CA and CB =", math.degrees(triangle1_angle_c))

print("\b")
#Triangle 2
# Angle AB and AC
triangle2_angle_a = math.acos((triangle2_ab.dot(triangle2_ac)) /
(triangle2_ab.norm() * triangle2_ac.norm()))
print("Triangle 2: Angle between AB and AC =", math.degrees(triangle2_angle_a))
# Angle BA and BC
triangle2_angle_b = math.acos((triangle2_ba.dot(triangle2_bc)) /
(triangle2_ba.norm() * triangle2_bc.norm()))
print("Triangle 2: Angle between BA and BC =", math.degrees(triangle2_angle_b))
# Angle CA and CB
triangle2_angle_c = math.acos((triangle2_ca.dot(triangle2_cb)) /
(triangle2_ca.norm() * triangle2_cb.norm()))
print("Triangle 2: Angle between CA and CB =", math.degrees(triangle2_angle_c))

print("\b")
#Triangle 3
# Angle AB and AC
triangle3_angle_a = math.acos((triangle3_ab.dot(triangle3_ac)) /
(triangle3_ab.norm() * triangle3_ac.norm()))
print("Triangle 3: Angle between AB and AC =", math.degrees(triangle3_angle_a))
# Angle BA and BC
triangle3_angle_b = math.acos((triangle3_ba.dot(triangle3_bc)) /
(triangle3_ba.norm() * triangle3_bc.norm()))
print("Triangle 3: Angle between BA and BC =", math.degrees(triangle3_angle_b))
# Angle CA and CB
triangle3_angle_c = math.acos((triangle3_ca.dot(triangle3_cb)) /
(triangle3_ca.norm() * triangle3_cb.norm()))
print("Triangle 3: Angle between CA and CB =", math.degrees(triangle3_angle_c))


print("\b")
#Triangle 4
# Angle AB and AC
triangle4_angle_a = math.acos((triangle4_ab.dot(triangle4_ac)) /
(triangle4_ab.norm() * triangle4_ac.norm()))
print("Triangle 4: Angle between AB and AC =", math.degrees(triangle4_angle_a))
# Angle BA and BC
triangle4_angle_b = math.acos((triangle4_ba.dot(triangle4_bc)) /
(triangle4_ba.norm() * triangle4_bc.norm()))
print("Triangle 4: Angle between BA and BC =", math.degrees(triangle4_angle_b))
# Angle CA and CB
triangle4_angle_c = math.acos((triangle4_ca.dot(triangle4_cb)) /
(triangle4_ca.norm() * triangle4_cb.norm()))
print("Triangle 4: Angle between CA and CB =", math.degrees(triangle4_angle_c))

print("c)")
#Triangle 1
triangle1_a_sp = vect.SphericalPolar(0, math.radians(0), math.radians(0))
triangle1_b_sp = vect.SphericalPolar(1, math.radians(0), math.radians(0))
triangle1_c_sp = vect.SphericalPolar(1, math.radians(90), math.radians(0))
triangle1_ab_sp = triangle1_b_sp - triangle1_a_sp
triangle1_ac_sp = triangle1_c_sp - triangle1_a_sp
triangle1_ba_sp = triangle1_a_sp - triangle1_b_sp
triangle1_bc_sp = triangle1_c_sp - triangle1_b_sp
triangle1_ca_sp = triangle1_a_sp - triangle1_c_sp
triangle1_cb_sp = triangle1_b_sp - triangle1_c_sp
# Area
triangle1_area_sp = 0.5 * (triangle1_ab_sp.cross(triangle1_ac_sp)).norm()
print("Area of Triangle 1 =", triangle1_area_sp)
# Angle AB and AC
triangle1_angle_a_sp = math.acos((triangle1_ab_sp.dot(triangle1_ac_sp)) /
(triangle1_ab_sp.norm() * triangle1_ac_sp.norm()))
print("Triangle 1: Angle between AB and AC =", math.degrees(triangle1_angle_a_sp))
# Angle BA and BC
triangle1_angle_b_sp = math.acos((triangle1_ba_sp.dot(triangle1_bc_sp)) /
(triangle1_ba_sp.norm() * triangle1_bc_sp.norm()))
print("Triangle 1: Angle between BA and BC =", math.degrees(triangle1_angle_b_sp))
# Angle CA and CB
triangle1_angle_c_sp = math.acos((triangle1_ca_sp.dot(triangle1_cb_sp)) /
(triangle1_ca_sp.norm() * triangle1_cb_sp.norm()))
print("Triangle 1: Angle between CA and CB =", math.degrees(triangle1_angle_c_sp))

#Triangle 2
triangle2_a_sp = vect.SphericalPolar(1, math.radians(0), math.radians(0))
triangle2_b_sp = vect.SphericalPolar(1, math.radians(90), math.radians(0))
triangle2_c_sp = vect.SphericalPolar(1, math.radians(90), math.radians(180))
triangle2_ab_sp = triangle2_b_sp - triangle2_a_sp
triangle2_ac_sp = triangle2_c_sp - triangle2_a_sp
triangle2_ba_sp = triangle2_a_sp - triangle2_b_sp
triangle2_bc_sp = triangle2_c_sp - triangle2_b_sp
triangle2_ca_sp = triangle2_a_sp - triangle2_c_sp
triangle2_cb_sp = triangle2_b_sp - triangle2_c_sp
triangle2_area_sp = 0.5 * (triangle2_ba_sp.cross(triangle2_bc_sp)).norm()
print("Area of Triangle 2 =", triangle2_area_sp)
# Angle AB and AC
triangle2_angle_a_sp = math.acos((triangle2_ab_sp.dot(triangle2_ac_sp)) /
(triangle2_ab_sp.norm() * triangle2_ac_sp.norm()))
print("Triangle 2: Angle between AB and AC =", math.degrees(triangle2_angle_a_sp))
# Angle BA and BC
triangle2_angle_b_sp = math.acos((triangle2_ba_sp.dot(triangle2_bc_sp)) /
(triangle2_ba_sp.norm() * triangle2_bc_sp.norm()))
print("Triangle 2: Angle between BA and BC =", math.degrees(triangle2_angle_b_sp))
# Angle CA and CB
triangle2_angle_c_sp = math.acos((triangle2_ca_sp.dot(triangle2_cb_sp)) /
(triangle2_ca_sp.norm() * triangle2_cb_sp.norm()))
print("Triangle 2: Angle between CA and CB =", math.degrees(triangle2_angle_c_sp))


#Triangle 3
triangle3_a_sp = vect.SphericalPolar(0, math.radians(0), math.radians(0))
triangle3_b_sp = vect.SphericalPolar(2, math.radians(0), math.radians(0))
triangle3_c_sp = vect.SphericalPolar(2, math.radians(90), math.radians(0))
triangle3_ab_sp = triangle3_b_sp - triangle3_a_sp
triangle3_ac_sp = triangle3_c_sp - triangle3_a_sp
triangle3_ba_sp = triangle3_a_sp - triangle3_b_sp
triangle3_bc_sp = triangle3_c_sp - triangle3_b_sp
triangle3_ca_sp = triangle3_a_sp - triangle3_c_sp
triangle3_cb_sp = triangle3_b_sp - triangle3_c_sp
# Area
triangle3_area_sp = 0.5 * (triangle3_ab_sp.cross(triangle3_ac_sp)).norm()
print("Area of Triangle 3 =", triangle3_area_sp)
# Angle AB and AC
triangle3_angle_a_sp = math.acos((triangle3_ab_sp.dot(triangle3_ac_sp)) /
(triangle3_ab_sp.norm() * triangle3_ac_sp.norm()))
print("Triangle 3: Angle between AB and AC =", math.degrees(triangle3_angle_a_sp))
# Angle BA and BC
triangle3_angle_b_sp = math.acos((triangle3_ba_sp.dot(triangle3_bc_sp)) /
(triangle3_ba_sp.norm() * triangle3_bc_sp.norm()))
print("Triangle 3: Angle between BA and BC =", math.degrees(triangle3_angle_b_sp))
# Angle CA and CB
triangle3_angle_c_sp = math.acos((triangle3_ca_sp.dot(triangle3_cb_sp)) /
(triangle3_ca_sp.norm() * triangle3_cb_sp.norm()))
print("Triangle 3: Angle between CA and CB =", math.degrees(triangle3_angle_c_sp))

#Triangle 4
triangle4_a_sp = vect.SphericalPolar(1, math.radians(90), math.radians(0))
triangle4_b_sp = vect.SphericalPolar(1, math.radians(90), math.radians(180))
triangle4_c_sp = vect.SphericalPolar(1, math.radians(90), math.radians(270))
triangle4_ab_sp = triangle4_b_sp - triangle4_a_sp
triangle4_ac_sp = triangle4_c_sp - triangle4_a_sp
triangle4_ba_sp = triangle4_a_sp - triangle4_b_sp
triangle4_bc_sp = triangle4_c_sp - triangle4_b_sp
triangle4_ca_sp = triangle4_a_sp - triangle4_c_sp
triangle4_cb_sp = triangle4_b_sp - triangle4_c_sp
# Area
triangle4_area_sp = 0.5 * (triangle4_ab_sp.cross(triangle4_ac_sp)).norm()
print("Area of Triangle 4 =", triangle4_area_sp)
# Angle AB and AC
triangle4_angle_a_sp = math.acos((triangle4_ab_sp.dot(triangle4_ac_sp)) /
(triangle4_ab_sp.norm() * triangle4_ac_sp.norm()))
print("Triangle 4: Angle between AB and AC =", math.degrees(triangle4_angle_a_sp))
# Angle BA and BC
triangle4_angle_b_sp = math.acos((triangle4_ba_sp.dot(triangle4_bc_sp)) /
(triangle4_ba_sp.norm() * triangle4_bc_sp.norm()))
print("Triangle 4: Angle between BA and BC =", math.degrees(triangle4_angle_b_sp))
# Angle CA and CB
triangle4_angle_c_sp = math.acos((triangle4_ca_sp.dot(triangle4_cb_sp)) /
(triangle4_ca_sp.norm() * triangle4_cb_sp.norm()))
print("Triangle 4: Angle between CA and CB =", math.degrees(triangle4_angle_c_sp))
