# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE
Title: post_processing.py
"""
# Basic packages:
from config_packages import np, math, plt, cm, data
from src import plotting
import src.quadrature as quadrature

def compute_gradient(U, mesh):
    gradient = np.zeros((mesh.ndof, 2))
    for i in range(len(U)):
        is_last_column = ((i + 1) % mesh.ndof_x == 0)
        is_last_row = (i >= (mesh.ndof_y - 1)*mesh.ndof_x)

        # dU/dx:
        if is_last_column:
            gradient[i][0] = 2*gradient[i-1][0] - gradient[i-2][0]
        else:
            gradient[i][0] = (U[i+1] - U[i])/mesh.h

        # dU/dy:
        if is_last_row:
            gradient[i][1] = 2*gradient[i - mesh.ndof_x][1] - gradient[i - 2*mesh.ndof_x][1]
        else:
            gradient[i][1] = (U[i + mesh.ndof_x] - U[i])/mesh.h

    return gradient

def interpolate_triangle(point, triangle_points, triangle_values):
    rhs = [point[0], point[1], 1]
    coordinates_matrix = [[triangle_points[0][0], triangle_points[1][0], triangle_points[2][0]], [triangle_points[0][1], triangle_points[1][1], triangle_points[2][1]], [1, 1, 1]]
    weights = np.linalg.solve(coordinates_matrix, rhs)

    interpolated_value = np.dot(weights, triangle_values)
    return interpolated_value

def error_u(point, triangle_points, triangle_values):
    return (data.u_analytical(point[0], point[1]) - interpolate_triangle(point, triangle_points, triangle_values))**2

def error_gradu(point, triangle_points, triangle_values):
    grad_dudx = [triangle_values[0][0], triangle_values[1][0], triangle_values[2][0]]
    grad_dudy = [triangle_values[0][1], triangle_values[1][1], triangle_values[2][1]]

    error_dudx = data.gradu_analytical(point[0], point[1])[0] - interpolate_triangle(point, triangle_points, grad_dudx)
    error_dudy = data.gradu_analytical(point[0], point[1])[1] - interpolate_triangle(point, triangle_points, grad_dudy)
    return error_dudx**2 + error_dudy**2

def compute_errors(U, mesh):
    # Analytical solution:
    U_analytical = np.zeros(mesh.ndof)
    gradU_analytical = np.zeros((mesh.ndof, 2))
    for index, node in enumerate(mesh.nodes):
        node_x, node_y = node[0], node[1]

        U_analytical[index] = data.u_analytical(node_x, node_y)
        gradU_analytical[index][0], gradU_analytical[index][1] = data.gradu_analytical(node_x, node_y)
    
    # Gradient computation of the numerical solution:
    gradU = compute_gradient(U, mesh)

    # Computation of L2, H1semi, H1 error norms:
    L2_error, H1semi_error = 0, 0
    for element in mesh.elements:
        det_J = np.linalg.det(element.J)
        points = [element.point_1, element.point_2, element.point_3]
        U_values = [U[element.nodes_indexes[0]], U[element.nodes_indexes[1]], U[element.nodes_indexes[2]]]
        gradU_values = [gradU[element.nodes_indexes[0]], gradU[element.nodes_indexes[1]], gradU[element.nodes_indexes[2]]]
        for quadrature_weight, quadrature_point in zip(quadrature.weights(), quadrature.points()):
            x_, y_ = quadrature.map(element, quadrature_point[0], quadrature_point[1])

            L2_error += abs(det_J)*quadrature_weight*error_u([x_, y_], points, U_values)
            H1semi_error += abs(det_J)*quadrature_weight*error_gradu([x_, y_], points, gradU_values)
    L2_error, H1semi_error = math.sqrt(L2_error), math.sqrt(H1semi_error)
    H1_error = math.sqrt(L2_error**2 + H1semi_error**2)

    plotting.plot_error(mesh, U, U_analytical, L2_error, H1_error)
    return L2_error, H1_error