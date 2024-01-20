# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE
Title: post_processing.py
"""
# Basic packages:
from config_packages import np, math, plt, cm, data

def x_map(epsilon, eta, element):
    x1 = element.x1
    x2 = element.x2
    x3 = element.x3
    return x1 + (x2 - x1)*epsilon + (x3 - x1)*eta

def y_map(epsilon, eta, element):
    y1 = element.y1
    y2 = element.y2
    y3 = element.y3
    return y1 + (y2 - y1)*epsilon + (y3 - y1)*eta

def compute_gradient(U, mesh):
    gradient_x = np.zeros(mesh.ndof)
    gradient_y = np.zeros(mesh.ndof)
    gradient = np.zeros((mesh.ndof, 2))
    # Calculation of du/dx:
    for i in range(len(U)):
        if (i + 1) % mesh.ndof_x == 0:
            # Extrapolation of the nodes on the right wall:
            gradient_x[i] = 2*gradient_x[i-1] - gradient_x[i-2]
        else:
            gradient_x[i] = (U[i+1] - U[i])/mesh.h

    # Calculation of du/dy:
    for i in range(len(U)):
        if i >= (mesh.ndof_y - 1)*mesh.ndof_x:
            # Extrapolation of the nodes on the top wall:
            gradient_y[i] = 2*gradient_y[i - mesh.ndof_x] - gradient_y[i - 2*mesh.ndof_x]
        else:
            gradient_y[i] = (U[i + mesh.ndof_x] - U[i])/mesh.h
    
    for i in range(len(U)):
        gradient[i][0] = gradient_x[i]
        gradient[i][1] = gradient_y[i]

    return gradient

def interpolate_triangle(point, triangle_points, triangle_values):
    rhs = [point[0], point[1], 1]

    coordinates_matrix = [[triangle_points[0][0], triangle_points[1][0], triangle_points[2][0]], [triangle_points[0][1], triangle_points[1][1], triangle_points[2][1]], [1, 1, 1]]
    interpolation_matrix = np.linalg.inv(coordinates_matrix)
    weights = np.dot(interpolation_matrix, rhs)

    interpolated_value = weights[0]*triangle_values[0] + weights[1]*triangle_values[1] + weights[2]*triangle_values[2]

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
    # Gaussian quadrature points:
    quad_points = [[0, 1/2], [1/2, 0], [1/2, 1/2]]
    quad_weights = [1/6, 1/6, 1/6]

    U_analytical = np.zeros(mesh.ndof)
    for i in range(mesh.ndof):
        U_analytical[i] = data.u_analytical(mesh.nodes[i][0], mesh.nodes[i][1])

    gradU_analytical = np.zeros((mesh.ndof, 2))
    for i in range(mesh.ndof):
        gradU_analytical[i][0] = data.gradu_analytical(mesh.nodes[i][0], mesh.nodes[i][1])[0]
        gradU_analytical[i][1] = data.gradu_analytical(mesh.nodes[i][0], mesh.nodes[i][1])[1]
    
    gradU = compute_gradient(U, mesh)

    L2_error = 0
    for element in mesh.elements:
        det_J = np.linalg.det(element.J)
        triangle_points = [element.point_1, element.point_2, element.point_3]
        triangle_values = [U[element.nodes_indexes[0]], U[element.nodes_indexes[1]], U[element.nodes_indexes[2]]]
        for quadrature_weight, quadrature_point in zip(quad_weights, quad_points):
            x_ = x_map(quadrature_point[0], quadrature_point[1], element)
            y_ = y_map(quadrature_point[0], quadrature_point[1], element)

            L2_error += abs(det_J)*quadrature_weight*error_u([x_, y_], triangle_points, triangle_values)
    L2_error = math.sqrt(L2_error)

    H1semi_error = 0
    for element in mesh.elements:
        det_J = np.linalg.det(element.J)
        triangle_points = [element.point_1, element.point_2, element.point_3]
        triangle_values = [gradU[element.nodes_indexes[0]], gradU[element.nodes_indexes[1]], gradU[element.nodes_indexes[2]]]
        for quadrature_weight, quadrature_point in zip(quad_weights, quad_points):
            x_ = x_map(quadrature_point[0], quadrature_point[1], element)
            y_ = y_map(quadrature_point[0], quadrature_point[1], element)

            H1semi_error += abs(det_J)*quadrature_weight*error_gradu([x_, y_], triangle_points, triangle_values)
    H1semi_error = math.sqrt(H1semi_error)

    H1_error = math.sqrt(L2_error**2 + H1semi_error**2)

    if data.plot_error == 'y':
       fig = plt.figure()
       ax = fig.add_subplot(projection='3d')
       ax.scatter(mesh.x_coord, mesh.y_coord, U, color ='black', marker = 'o', label = "Finite Element Solution")
       ax.scatter(mesh.x_coord, mesh.y_coord, U_analytical, color ='purple', marker = 'x', label = "Analytical Solution")
       ax.legend()
       plt.title('Finite Element Solver')
       ax.set_xlabel('x')
       ax.set_ylabel('y')
       ax.set_zlabel('u(x,y)')
       plt.figtext(0.15, 0.83, "L2 Error = " + '%.6f' %L2_error)
       plt.figtext(0.15, 0.75, "H1 Error = " + '%.6f' %H1_error)
       plt.show()

    return L2_error, H1_error