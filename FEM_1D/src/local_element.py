# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE FEM 1D
Title: local_element.py
"""
# Basic packages:
from config_packages import np, math, plt, cm, data

# Gaussian quadrature points:
quad_points = [-math.sqrt(3/5), 0 , + math.sqrt(3/5)]
quad_weights = [5/9, 8/9 , 5/9]

class Local_element:
    def __init__(self, point_1, point_2):
        self.x1 = point_1
        self.x2 = point_2
        self.T = [[1, point_1[0]],[1, point_2[0]]]

    def phi(self, n, x_):
        # Coefficients of the basis function:
        a0 = np.linalg.inv(self.T)[0][n]
        a1 = np.linalg.inv(self.T)[1][n]

        return a0 + a1*x_

    def dphi(self, n, x_):
        # Coefficients of the basis function:
        a0 = np.linalg.inv(self.T)[0][n]
        a1 = np.linalg.inv(self.T)[1][n]

        return a1

    def x_map(self, epsilon):
        return 0.5*(self.x1 + self.x2) + 0.5*(self.x2 - self.x1)*epsilon

    def A_local(self):
        A_local_ = np.zeros([2, 2])
        K_local = np.zeros([2, 2]) # Local stiffness matrix
        V_local = np.zeros([2, 2]) # Local advection matrix
        M_local = np.zeros([2, 2]) # Local mass matrix
        det_J = 0.5*(self.x2 - self.x1)
        for m in range(2):
            for n in range(2):
                for quadrature_weight, quadrature_point in zip(quad_weights, quad_points):
                    x_ = self.x_map(quadrature_point)
                    phi_m = self.phi(m, x_)
                    phi_n = self.phi(n, x_)
                    dphi_m = self.dphi(m, x_)
                    dphi_n = self.dphi(n, x_)

                    K_local[m][n] += abs(det_J)*quadrature_weight*data.mu(x_)*dphi_m*dphi_n
                    V_local[m][n] += abs(det_J)*quadrature_weight*data.beta(x_)*dphi_n*phi_m
                    M_local[m][n] += abs(det_J)*quadrature_weight*data.sigma(x_)*phi_m*phi_n

                    A_local_ = K_local + V_local + M_local
        return A_local_

    def F_local(self):
        F_local_ = np.zeros(2)
        det_J = 0.5*(self.x2 - self.x1)
        for m in range(2):
            for quadrature_weight, quadrature_point in zip(quad_weights, quad_points):
                x_ = self.x_map(quadrature_point)
                phi_m = self.phi(m, x_)

                F_local_[m] += abs(det_J)*quadrature_weight*data.f(x_)*phi_m
        return F_local_