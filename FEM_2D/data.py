# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE
Title: data.py
"""

import numpy as np 
import math

# Plot options:
plot_solution = 'n'
plot_error = 'y'

# Parameters:
Lx = 10**5
Ly = 2*math.pi*10**4
H = 200
W = 0.3*10**(-7)
R = 0.6*10**(-3)
beta = 5*10**(-10)
#beta = 0

# Rectangular domain definition:
x1 = 0
y1 = 0
x2 = Lx
y2 = Ly
h = Lx/10

alpha = H*beta/R
gamma = W*math.pi/(R*Ly)

def f(x, y):
    return gamma*math.sin(math.pi*y/Ly)

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
    b0 = - alpha
    b1 = 0
    return [b0, b1]

def sigma(x,y):
    return 0

def u_an(x,y):
    lambda_1 = -alpha*0.5 + math.sqrt((alpha*0.5)**2 + (math.pi/Ly)**2)
    lambda_2 = -alpha*0.5 - math.sqrt((alpha*0.5)**2 + (math.pi/Ly)**2) 
    p = (1 - math.exp(lambda_2*Lx))/(math.exp(lambda_1*Lx) - math.exp(lambda_2*Lx))

    psi_an = - gamma*((Ly/math.pi)**2)*(p*math.exp(lambda_1*x) + (1 - p)*math.exp(lambda_2*x) - 1)*math.sin(math.pi*y/Ly)
    return psi_an