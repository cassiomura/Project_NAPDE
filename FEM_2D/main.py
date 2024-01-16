# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE
Title: main.py
"""
# Basic packages:
from config_packages import np, math, plt, cm, data, module_name

# Custom packages:
from src import mesh_generation, boundary_conditions, post_processing

#0. Read data:
print("============================================================")
print("(0/6) Reading data from " + module_name +  ".py ...")
print("============================================================")

#1. Mesh generation:
print("============================================================")
print("(1/6) Generating mesh ...")
print("============================================================")
mesh = mesh_generation.Mesh([data.x1, data.x2], [data.y1, data.y2], data.h)

#2. Assemble of global matrices and right hand side:
print("============================================================")
print("(2/6) Assembling global matrices and right hand side ...")
print("============================================================")
A = np.zeros([mesh.ndof, mesh.ndof])
element = 0 
for local_element in mesh.elements:
    for m in range(3):
        for n in range(3):
            A[mesh.elements_nodes[element][m]][mesh.elements_nodes[element][n]] += local_element.A_local()[m][n] + local_element.V_local()[m][n] + local_element.M_local()[m][n]
    element = element + 1 

F = np.zeros(mesh.ndof)
element = 0
for local_element in mesh.elements:
    for m in range(3):
        F[mesh.elements_nodes[element][m]] += local_element.F_local()[m]
    element = element + 1

# 3. Impose boundary conditions:
print("============================================================")
print("(3/6) Imposing boundary conditions ...")
print("============================================================")
A, F, g = boundary_conditions.impose_boundary_conditions(A, F, mesh)

# 4. Solve the algebraic problem:
print("============================================================")
print("(4/6) Solving the algebric problem ...")
print("============================================================")
U = np.linalg.solve(A, F)
# Lifting operation:
U = U + g

# 5. Plotting the solution:
print("============================================================")
print("(5/6) Plotting the solution ...")
print("============================================================")
# Plot the solution
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
print("(6/6) Computing errors ...")
print("============================================================")
err = post_processing.compute_errors(U, mesh)

# 7. Finish the program:
print("============================================================")
print("Program completed successfully! ")
print("============================================================")
exit()