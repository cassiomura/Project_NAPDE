# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE
Title: main.py
"""

# Basic packages:
from config_packages import np, math, plt, cm

# Custom packages:
import mesh_generation
import boundary_conditions
import post_processing
import data

#1. Mesh generation:
print("============================================================")
print(" Generating mesh ...")
print("============================================================")
mesh = mesh_generation.Mesh([data.x1, data.x2], [data.y1, data.y2], data.h)

#2. Assemble of global matrices and right hand side:
print("============================================================")
print(" Assembling global matrices and right hand side ...")
print("============================================================")
A = np.zeros([mesh.ndof, mesh.ndof])
elem_index = 0 
for local_element in mesh.elements:
    for m in range(3):
        for n in range(3):
            A[mesh.elements_nodes_indexes[elem_index][m]][mesh.elements_nodes_indexes[elem_index][n]] += local_element.A_local()[m][n] + local_element.V_local()[m][n] + local_element.M_local()[m][n]
    elem_index = elem_index + 1 

F = np.zeros(mesh.ndof)
elem_index = 0
for local_element in mesh.elements:
    for m in range(3):
        F[mesh.elements_nodes_indexes[elem_index][m]] += local_element.F_local()[m]
    elem_index = elem_index + 1

# 3. Impose boundary conditions:
print("============================================================")
print(" Imposing boundary conditions ...")
print("============================================================")
A, F, g = boundary_conditions.impose_boundary_conditions(A, F, mesh)

# 4. Solve the algebraic problem:
print("============================================================")
print(" Solving the algebric problem ...")
print("============================================================")
U = np.linalg.solve(A, F)
# Lifting:
U = U + g

# 5. Post Processing:
print("============================================================")
print(" Post-processing the solution ...")
print("============================================================")
if data.plot_solution == 'y':
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    surf = ax.plot_trisurf(mesh.x_coord, mesh.y_coord, U, cmap=cm.coolwarm)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('u(x,y)')
    plt.title('Finite Element Solver')
    plt.show()

# 6. Computing the error:
print("============================================================")
print(" Computing errors ...")
print("============================================================")
err = post_processing.compute_errors(U, mesh)