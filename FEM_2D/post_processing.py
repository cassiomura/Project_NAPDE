# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE
Title: Post-Processing.py
"""
import numpy as np 
import math
import matplotlib.pyplot as plt
from matplotlib import cm

from data import u_an, plot_error

def compute_errors(U, mesh):
    U_an = np.zeros(mesh.ndof)
    for i in range(mesh.ndof):
        U_an[i] = u_an(mesh.nodes[i][0], mesh.nodes[i][1])

    err = 0
    for i in range(mesh.ndof):
        err += (1/mesh.ndof)*(U[i] - U_an[i])**2
    err = math.sqrt(err)

    if plot_error == 'y':
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