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

def interpolate_triangle(point, triangle_points, triangle_values):
    rhs = [point[0], point[1], 1]

    coordinates_matrix = [[triangle_points[0][0], triangle_points[1][0], triangle_points[2][0]], [triangle_points[0][1], triangle_points[1][1], triangle_points[2][1]], [1, 1, 1]]
    interpolation_matrix = np.linalg.inv(coordinates_matrix)
    weights = np.dot(interpolation_matrix, rhs)

    interpolated_value = weights[0]*triangle_values[0] + weights[1]*triangle_values[1] + weights[2]*triangle_values[2]

    return interpolated_value

def error_u(point, triangle_points, triangle_values):
    return (data.u_analytical(point[0], point[1]) - interpolate_triangle(point, triangle_points, triangle_values))**2

def compute_errors(U, mesh):
    # Gaussian quadrature points:
    quad_points = [[0, 1/2], [1/2, 0], [1/2, 1/2]]
    quad_weights = [1/6, 1/6, 1/6]

    U_analytical = np.zeros(mesh.ndof)
    for i in range(mesh.ndof):
        U_analytical[i] = data.u_analytical(mesh.nodes[i][0], mesh.nodes[i][1])
    
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
       plt.show()

    return L2_error


#def compute_errors(U, mesh):
#    U_an = np.zeros(mesh.ndof)
#    for i in range(mesh.ndof):
#        U_an[i] = data.u_an(mesh.nodes[i][0], mesh.nodes[i][1])

#    err = 0
#    for i in range(mesh.ndof):
#        err += (1/mesh.ndof)*(U[i] - U_an[i])**2
#    err = math.sqrt(err)

#    if data.plot_error == 'y':
#        fig = plt.figure()
#        ax = fig.add_subplot(projection='3d')
#        ax.scatter(mesh.x_coord, mesh.y_coord, U, color ='black', marker = 'o', label = "Finite Element Solution")
#        ax.scatter(mesh.x_coord, mesh.y_coord, U_an, color ='purple', marker = 'x', label = "Analytical Solution")
#        ax.legend()
#        plt.title('Finite Element Solver')
#        ax.set_xlabel('x')
#        ax.set_ylabel('y')
#        ax.set_zlabel('u(x,y)')
#        plt.figtext(0.15, 0.83, "Error = " + '%.6f' %err)
#        plt.show()

#    return err