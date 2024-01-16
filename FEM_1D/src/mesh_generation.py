# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE FEM 1D
Title: mesh_generation.py
"""
# Basic packages:
from config_packages import np, math
from src.local_element import Local_element

class Mesh:
    def __init__(self, x_domain, h):
        self.h = h

        ndof_ = (math.floor((x_domain[1] - x_domain[0])/h) + 1)
        self.ndof = ndof_
        
        nodes_ = np.zeros((ndof_, 1))
        dof_counter = 0
        for i in range(ndof_):
            nodes_[dof_counter][0] = i*h + x_domain[0] 
            dof_counter += 1
        self.nodes = nodes_

        coord_ = np.zeros(len(nodes_))
        for i in range(len(nodes_)):
            coord_[i] = nodes_[i]
        self.coord = coord_

        elements_nodes_ = []
        for i in range(ndof_ - 1):
            elements_nodes_.append([i, i + 1])
        self.elements_nodes = elements_nodes_

        elements_ = []
        for element_indexes in elements_nodes_:
            node1 = nodes_[element_indexes[0]]
            node2 = nodes_[element_indexes[1]]
            elements_.append(Local_element(node1, node2))
        self.elements = elements_