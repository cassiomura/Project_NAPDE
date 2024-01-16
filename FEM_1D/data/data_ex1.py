# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE
Title: data_ex1.py
"""
import numpy as np 
import math

#https://www.wolframalpha.com/input?i2d=true&i=-+D%5Bu%2C%7Bx%2C2%7D%5D+%2B+u+%3D+sin%5C%2840%29x%5C%2841%29%5C%2844%29+u%5C%2840%290%5C%2841%29%3D0%5C%2844%29+u%5C%2840%292+pi%5C%2841%29+%3D+0

# Domain definition:
x1 = 0
x2 = 2*math.pi
h = 0.01

# Plot options:
plot_solution = 'n'
plot_error = 'y'

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
    return 0.5*math.sin(x)