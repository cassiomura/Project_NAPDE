# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE
Title: mesh_generation.py
"""
# Basic packages:
from config_packages import np, math
from src.local_element import Local_element

class Mesh:
    def __init__(self, x_domain, y_domain, h):
        self.h = h

        ndof_x_ = (math.floor((x_domain[1] - x_domain[0])/h) + 1)
        self.ndof_x = ndof_x_

        ndof_y_ = (math.floor((y_domain[1] - y_domain[0])/h) + 1)
        self.ndof_y = ndof_y_

        ndof_ = ndof_x_*ndof_y_
        self.ndof = ndof_
        
        nodes_ = np.zeros((ndof_, 2))
        dof_counter = 0
        for i in range(ndof_y_):
            for j in range(ndof_x_):
                nodes_[dof_counter][0] = j*h + x_domain[0] 
                nodes_[dof_counter][1] = i*h + y_domain[0]
                dof_counter += 1
        self.nodes = nodes_

        x_coord_ = np.zeros(len(nodes_))
        y_coord_ = np.zeros(len(nodes_))
        for i in range(len(nodes_)):
            x_coord_[i] = nodes_[i][0]
            y_coord_[i] = nodes_[i][1]
        self.x_coord = x_coord_
        self.y_coord = y_coord_

        elements_nodes_ = []
        for j in range(ndof_y_ - 1):
            for i in range(ndof_x_ - 1):
                p = ndof_x_*j + i
                elements_nodes_.append([p, p + ndof_x_ + 1, p + ndof_x_])
                elements_nodes_.append([p, p + ndof_x_ + 1, p + 1])
        self.elements_nodes = elements_nodes_

        elements_ = []
        for element_indexes in elements_nodes_:
            node1 = nodes_[element_indexes[0]]
            node2 = nodes_[element_indexes[1]]
            node3 = nodes_[element_indexes[2]]
            elements_.append(Local_element(node1, node2, node3, element_indexes))
        self.elements = elements_