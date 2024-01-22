# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE
Title: post_processing.py
"""
# Basic packages:
from config_packages import np, math, plt, cm, data

def plot_solution(mesh, U):
    if data.plot_solution == 'y':
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        surf = ax.plot_trisurf(mesh.x_coord, mesh.y_coord, U, cmap=cm.coolwarm)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('u(x,y)')
        plt.title('Finite Element Solver')
        plt.draw()

def plot_error(mesh, U, U_analytical, L2_error, H1_error):
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