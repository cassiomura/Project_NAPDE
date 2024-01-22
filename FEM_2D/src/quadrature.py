# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE
Title: quadrature.py
"""
# Basic packages:
from config_packages import np, math, plt, cm, data

def weights():
    return [1/6, 1/6, 1/6]

def points():
    return [[0, 1/2], [1/2, 0], [1/2, 1/2]]

def map(element, epsilon, eta):
    x1, y1 = element.point_1[0], element.point_1[1]
    x2, y2 = element.point_2[0], element.point_2[1]
    x3, y3 = element.point_3[0], element.point_3[1]

    x_mapped = x1 + (x2 - x1)*epsilon + (x3 - x1)*eta
    y_mapped = y1 + (y2 - y1)*epsilon + (y3 - y1)*eta
    return x_mapped, y_mapped