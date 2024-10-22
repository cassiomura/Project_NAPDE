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
y1 = 0
x2 = 1
y2 = 1
h = 0.05

# Plot options:
plot_solution = 'y'
plot_error = 'y'

def f(x, y):
    return 8*math.pi**2*math.sin(2*math.pi*x)*math.sin(2*math.pi*y)

def g1(x,y):
    return 0

def g2(x,y):
    return 0

def g3(x,y):
    return 0

def g4(x,y):
    return 0

def mu(x, y):
    return 1

def beta(x, y):
    b0 = 0
    b1 = 0
    return [b0, b1]

def sigma(x,y):
    return 0

def u_analytical(x,y):
    return math.sin(2*math.pi*x)*math.sin(2*math.pi*y)

def gradu_analytical(x,y):
    dudx = 2*math.pi*math.cos(2*math.pi*x)*math.sin(2*math.pi*y)
    dudy = 2*math.pi*math.sin(2*math.pi*x)*math.cos(2*math.pi*y)
    return [dudx, dudy]