# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE
Title: data_ex1.py
"""
import numpy as np 
import math

# Domain definition:
x1 = 0
x2 = 2*math.pi
h = 0.1

# Plot options:
plot_solution = 'y'
plot_error = 'n'

def f(x):
    return math.sin(x)

def g1(x):
    return 0

def g2(x):
    return 0


def mu(x):
    return 1

def beta(x):
    return 0

def sigma(x):
    return 1

def u_an(x):
    return math.sin(x)