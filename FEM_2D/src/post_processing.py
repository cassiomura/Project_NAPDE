# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE
Title: post_processing.py
"""
# Basic packages:
from config_packages import np, math, plt, cm, data

def compute_errors(U, mesh):
    U_an = np.zeros(mesh.ndof)
    for i in range(mesh.ndof):
        U_an[i] = data.u_an(mesh.nodes[i][0], mesh.nodes[i][1])

    err = 0
    for i in range(mesh.ndof):
        err += (1/mesh.ndof)*(U[i] - U_an[i])**2
    err = math.sqrt(err)

    if data.plot_error == 'y':
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(mesh.x_coord, mesh.y_coord, U, color ='black', marker = 'o', label = "Finite Element Solution")
        ax.scatter(mesh.x_coord, mesh.y_coord, U_an, color ='purple', marker = 'x', label = "Analytical Solution")
        ax.legend()
        plt.title('Finite Element Solver')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('u(x,y)')
        plt.figtext(0.15, 0.83, "Error = " + '%.6f' %err)
        plt.show()

    return err