# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE FEM 1D
Title: post_processing.py
"""
# Basic packages:
from config_packages import np, math, plt, cm, data

# Should have been called interpolate:
def u_h(x , x1, y1, x2, y2):
        return (x - x1)*(y2 - y1)/(x2 - x1) + y1

def err(x, x1, y1, x2, y2):
    return (data.u_analytical(x) - u_h(x, x1, y1, x2, y2))**2

def err_grad(x, x1, y1, x2, y2):
    return (data.grad_u_analytical(x) - u_h(x, x1, y1, x2, y2))**2

def compute_gradient(U, mesh):
    gradient_U = np.zeros(mesh.ndof)
    for i in range(len(U) - 1):
        gradient_U[i] = (U[i+1] - U[i])/mesh.h
    # Extrapolation for computing the last element:
    gradient_U[-1] = 2*gradient_U[-2] - gradient_U[-3] 

    return gradient_U

def compute_errors(U, mesh):
    # Gaussian quadrature points:
    quad_points = [-math.sqrt(3/5), 0 , + math.sqrt(3/5)]
    quad_weights = [5/9, 8/9 , 5/9]

    # Computation of the gradient:
    grad_U = compute_gradient(U, mesh)

    U_analytical = np.zeros(mesh.ndof)
    grad_U_analytical = np.zeros(mesh.ndof)
    for i in range(mesh.ndof):
        U_analytical[i] = data.u_analytical(mesh.nodes[i])
        grad_U_analytical[i] = data.grad_u_analytical(mesh.nodes[i])

    L2_error = 0
    for element in mesh.elements:
        for quadrature_weight, quadrature_point in zip(quad_weights, quad_points):
            x_map = 0.5*(element.x2 - element.x1)*quadrature_point + 0.5*(element.x2 + element.x1)
            L2_error += 0.5*(element.x2 - element.x1)*quadrature_weight*err(x_map, element.x1, U[element.index_node_1], element.x2, U[element.index_node_2])
    L2_error = math.sqrt(L2_error)

    H1_semi_error = 0
    for element in mesh.elements:
        for quadrature_weight, quadrature_point in zip(quad_weights, quad_points):
            x_map = 0.5*(element.x2 - element.x1)*quadrature_point + 0.5*(element.x2 + element.x1)
            H1_semi_error += 0.5*(element.x2 - element.x1)*quadrature_weight*err_grad(x_map, element.x1, grad_U[element.index_node_1], element.x2, grad_U[element.index_node_2])
    H1_semi_error = math.sqrt(H1_semi_error)

    H1_error = math.sqrt(L2_error**2 + H1_semi_error**2)

    if data.plot_error == 'y':
        fig = plt.figure()
        plt.scatter(mesh.coord, U, color ='black', marker = 'o', label = "Finite Element Solution")
        plt.scatter(mesh.coord, U_analytical, color ='purple', marker = 'x', label = "Analytical Solution")
        plt.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1)
        plt.title('Finite Element Solver')
        plt.xlabel('x')
        plt.ylabel('u(x)')
        plt.figtext(0.15, 0.83, "L2-error = " + '%.6f' %L2_error)
        plt.show()
        
    return L2_error, H1_error