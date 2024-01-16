# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE FEM 1D
Title: post_processing.py
"""
# Basic packages:
from config_packages import np, math, plt, cm, data

def u_h(x , x1, y1, x2, y2):
        return (x - x1)*(y2 - y1)/(x2 - x1) + y1

def err(x, x1, y1, x2, y2):
        return (data.u_an(x) - u_h(x, x1, y1, x2, y2))**2

def compute_errors(U, mesh):
    # Gaussian quadrature points:
    quad_points = [-math.sqrt(3/5), 0 , + math.sqrt(3/5)]
    quad_weights = [5/9, 8/9 , 5/9]

    U_an = np.zeros(mesh.ndof)
    for i in range(mesh.ndof):
        U_an[i] = data.u_an(mesh.nodes[i])

    L2_error = 0
    for element in mesh.elements:
        for quadrature_weight, quadrature_point in zip(quad_weights, quad_points):
            x_map = 0.5*(element.x2 - element.x1)*quadrature_point + 0.5*(element.x2 + element.x1)
            L2_error += 0.5*(element.x2 - element.x1)*quadrature_weight*err(x_map, element.x1, U[element.index_node_1], element.x2, U[element.index_node_2])
    L2_error = math.sqrt(L2_error)
    #err = 0
    #for i in range(mesh.ndof):
    #    err += (1/mesh.ndof)*(U[i] - U_an[i])**2
    #err = math.sqrt(err)

    #if data.plot_error == 'y':
    #    fig = plt.figure()
    #    ax = fig.add_subplot(projection='3d')
    #    ax.scatter(mesh.x_coord, mesh.y_coord, U, color ='black', marker = 'o', label = "Finite Element Solution")
    #    ax.scatter(mesh.x_coord, mesh.y_coord, U_an, color ='purple', marker = 'x', label = "Analytical Solution")
    #    ax.legend()
    #    plt.title('Finite Element Solver')
    #    ax.set_xlabel('x')
    #    ax.set_ylabel('y')
    #    ax.set_zlabel('u(x,y)')
    #    plt.figtext(0.15, 0.83, "Error = " + '%.6f' %err)
    #    plt.show()

    return L2_error