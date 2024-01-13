# Author: Cássio Murakami
import numpy as np 
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from data import f, mu, beta, sigma

# Gaussian quadrature points:
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
        det_J = np.linalg.det(self.J)
        for m in range(3):
            for n in range(3):
                for q in range(len(quad_points)):
                    x_ = self.x_map(quad_points[q][0], quad_points[q][1])
                    y_ = self.y_map(quad_points[q][0], quad_points[q][1])
                    dphi_m = self.dphi(m, x_, y_)
                    dphi_n = self.dphi(n, x_, y_)

                    A_local_[m][n] += abs(det_J)*quad_weights[q]*mu(x_, y_)*np.dot(dphi_m, dphi_n)
        return A_local_

    def V_local(self):
        V_local_ = np.zeros([3, 3])
        det_J = np.linalg.det(self.J)
        for m in range(3):
            for n in range(3):
                for q in range(len(quad_points)):
                    x_ = self.x_map(quad_points[q][0], quad_points[q][1])
                    y_ = self.y_map(quad_points[q][0], quad_points[q][1])
                    phi_m = self.phi(m, x_, y_)
                    dphi_n = self.dphi(n, x_, y_)

                    V_local_[m][n] += abs(det_J)*quad_weights[q]*np.dot(beta(x_, y_), dphi_n)*phi_m
        return V_local_

    def M_local(self):
        M_local_ = np.zeros([3, 3])
        det_J = np.linalg.det(self.J)
        for m in range(3):
            for n in range(3):
                for q in range(len(quad_points)):
                    x_ = self.x_map(quad_points[q][0], quad_points[q][1])
                    y_ = self.y_map(quad_points[q][0], quad_points[q][1])
                    phi_m = self.phi(m, x_, y_)
                    phi_n = self.phi(n, x_, y_)

                    M_local_[m][n] += abs(det_J)*quad_weights[q]*sigma(x_, y_)*phi_m*phi_n
        return M_local_

    def F_local(self):
        F_local_ = np.zeros(3)
        det_J = np.linalg.det(self.J)
        for m in range(3):
            for q in range(len(quad_points)):
                x_ = self.x_map(quad_points[q][0], quad_points[q][1])
                y_ = self.y_map(quad_points[q][0], quad_points[q][1])
                phi_m = self.phi(m, x_, y_)

                F_local_[m] += abs(det_J)*quad_weights[q]*f(x_, y_)*phi_m
        return F_local_