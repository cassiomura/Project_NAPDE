# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE
Title: Boundary-condition.py
"""
# Basic packages:
from config_packages import np, math, plt, cm, data

def impose_boundary_conditions(A, F, mesh):
    g = np.zeros(mesh.ndof)
    A, F, g = impose_Dirichlet("left", A, F, g, mesh)
    A, F, g = impose_Dirichlet("bottom", A, F, g, mesh)
    A, F, g = impose_Dirichlet("right", A, F, g, mesh)
    A, F, g = impose_Dirichlet("top", A, F, g, mesh)
    return A, F, g

def impose_Dirichlet(wall, A_, F_, g_, mesh_):
    bc_index_list = []
    g = g_
    g_wall = np.zeros(mesh_.ndof)
    for i in range(len(mesh_.nodes)):
        if (wall == "left" and mesh_.nodes[i][0] == mesh_.nodes[0][0]):
            bc_index_list.append(i)
            g_wall[i] = data.g1(mesh_.nodes[i][0], mesh_.nodes[i][1])
            g[i] = g_wall[i]
        elif (wall == "right" and mesh_.nodes[i][0] == mesh_.nodes[-1][0]):
            bc_index_list.append(i)
            g_wall[i] = data.g2(mesh_.nodes[i][0], mesh_.nodes[i][1])
            g[i] = g_wall[i]
        elif (wall == "bottom" and mesh_.nodes[i][1] == mesh_.nodes[0][1]):
            bc_index_list.append(i)
            g_wall[i] = data.g3(mesh_.nodes[i][0], mesh_.nodes[i][1])
            g[i] = g_wall[i]
        elif (wall == "top" and mesh_.nodes[i][1] == mesh_.nodes[-1][1]):
            bc_index_list.append(i)
            g_wall[i] = data.g4(mesh_.nodes[i][0], mesh_.nodes[i][1])
            g[i] = g_wall[i]

    # Lifting:
    F = F_ - A_.dot(g_wall)

    A = A_
    for bc_index in bc_index_list:
        for i in range(len(A_)):
            A[bc_index][i] = 0
            A[i][bc_index] = 0
        A[bc_index][bc_index] = 1
        F[bc_index] = 0 
    
    return A, F, g