# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE FEM_1D
Title: local_element.py
"""
# Basic packages:
from config_packages import np, math, plt, cm, data

# Gaussian quadrature points:
#zeta_list = [-math.sqrt(3/5), 0 , + math.sqrt(3/5)]
#weight_list = [5/9, 8/9 , 5/9]

quad_points = [[0, 1/2], [1/2, 0], [1/2, 1/2]]
quad_weights = [1/6, 1/6, 1/6]

class Local_element:
    def __init__(self, point_1, point_2, point_3):
        self.x1 = point_1[0]
        self.y1 = point_1[1]
        self.x2 = point_2[0]
        self.y2 = point_2[1]
        self.x3 = point_3[0]
        self.y3 = point_3[1]
        self.J = [point_2[0] - point_1[0], point_2[1]- point_1[1]],[point_3[0] - point_1[0], point_3[1] -point_1[1]]
        self.T = [[1, point_1[0], point_1[1]],[1, point_2[0], point_2[1]], [1, point_3[0], point_3[1]]]

    def phi(self, n, x_, y_):
        # Coefficients of the basis function:
        a0 = np.linalg.inv(self.T)[0][n]
        a1 = np.linalg.inv(self.T)[1][n]
        a2 = np.linalg.inv(self.T)[2][n]

        return a0 + a1*x_ + a2*y_

    def dphi(self, n, x_, y_):
        # Coefficients of the basis function:
        a0 = np.linalg.inv(self.T)[0][n]
        a1 = np.linalg.inv(self.T)[1][n]
        a2 = np.linalg.inv(self.T)[2][n]

        return [a1, a2]

    def x_map(self, epsilon, eta):
        return self.x1 + (self.x2 - self.x1)*epsilon + (self.x3 - self.x1)*eta

    def y_map(self, epsilon, eta):
        return self.y1 + (self.y2 - self.y1)*epsilon + (self.y3 - self.y1)*eta

    def A_local(self):
        A_local_ = np.zeros([3, 3])

        K_local = np.zeros([3, 3]) # Local stiffness matrix
        V_local = np.zeros([3, 3]) # Local advection matrix
        M_local = np.zeros([3, 3]) # Local mass matrix

        det_J = np.linalg.det(self.J)
        for m in range(3):
            for n in range(3):
                for quadrature_weight, quadrature_point in zip(quad_weights, quad_points):
                    x_ = self.x_map(quadrature_point[0], quadrature_point[1])
                    y_ = self.y_map(quadrature_point[0], quadrature_point[1])
                    phi_m = self.phi(m, x_, y_)
                    phi_n = self.phi(n, x_, y_)
                    dphi_m = self.dphi(m, x_, y_)
                    dphi_n = self.dphi(n, x_, y_)

                    K_local[m][n] += abs(det_J)*quadrature_weight*data.mu(x_, y_)*np.dot(dphi_m, dphi_n)
                    V_local[m][n] += abs(det_J)*quadrature_weight*np.dot(data.beta(x_, y_), dphi_n)*phi_m
                    M_local[m][n] += abs(det_J)*quadrature_weight*data.sigma(x_, y_)*phi_m*phi_n

                    A_local_ = K_local + V_local + M_local
        return A_local_

    def F_local(self):
        F_local_ = np.zeros(3)
        det_J = np.linalg.det(self.J)
        for m in range(3):
            for quadrature_weight, quadrature_point in zip(quad_weights, quad_points):
                x_ = self.x_map(quadrature_point[0], quadrature_point[1])
                y_ = self.y_map(quadrature_point[0], quadrature_point[1])
                phi_m = self.phi(m, x_, y_)

                F_local_[m] += abs(det_J)*quadrature_weight*data.f(x_, y_)*phi_m
        return F_local_