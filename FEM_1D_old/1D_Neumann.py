# Author: Cássio Murakami
# Finite Element Solver: 1D Laplacian
import numpy as np 
import math
import matplotlib.pyplot as plt

def f(x):
    return math.sin(x)

def sigma(x):
    return 1

def impose_Neumann(F, index, h):
    F_ = F
    F_[index] = F_[index] + h
    return F_

def impose_Dirichlet(A, F, index, g):
    A_ = A
    F_ = F
    for i in range(len(A)):
        A_[index][i] = 0
    A_[index][index] = 1
    F_[index] = g
    return (A_, F_)

def phi_i(x, a, b):
    return (b - x)/(b - a)

def dphi_i(x, a, b):
    return -1/(b - a)

def phi_j(x, a, b):
    return (x - a)/(b - a)

def dphi_j(x, a, b):
    return 1/(b - a)

def F_i(x, a, b):
    return phi_i(x, a, b)*f(x)

def F_j(x, a, b):
    return phi_j(x, a, b)*f(x)

def A_ii(x, a, b):
    return dphi_i(x, a, b)*dphi_i(x, a, b)

def A_ij(x, a, b):
    return dphi_i(x, a, b)*dphi_j(x, a, b)

def A_jj(x, a, b):
    return dphi_j(x, a, b)*dphi_j(x, a, b)

def M_ii(x, a, b):
    return sigma(x)*phi_i(x, a, b)*phi_i(x, a, b)

def M_ij(x, a, b):
    return sigma(x)*phi_i(x, a, b)*phi_j(x, a, b)

def M_jj(x, a, b):
    return sigma(x)*phi_j(x, a, b)*phi_j(x, a, b)

def uh_interpolate(x, x1, x2, u1, u2):
    return (x - x1)(u2 - u1)/(x2 - x1) + u1

def gauss_integration(func, a, b):
    # Gaussian quadrature parameters:
    zeta_list = [-math.sqrt(3/5), 0 , + math.sqrt(3/5)]
    weight_list = [5/9, 8/9 , 5/9]

    # Gaussian integration loop:
    result = 0
    for q in range(len(zeta_list)):
        x_ = zeta_list[q]*(b - a)/2 + (a + b)/2
        result = result + 0.5*(b-a)*weight_list[q]*func(x_, a, b)
    return result

#0. Problem paramenters:
# Domain boundaries
x1 = 0
x2 = 2*math.pi
# Dirichlet boundary conditions:
g1 = 0
g2 = 0
# Mesh size
h = 0.01

#1. Mesh generation:
dof = math.floor((x2 - x1)/h) + 1
x = np.zeros(dof)
for i in range(dof):
    x[i] = i*h + x1

class Local_element:
    def __init__(self, index_i, index_j):
        self.index_i = index_i
        self.index_j = index_j

    def A_local(self):
        A_local_ = np.zeros([2, 2])
        A_local_[0][0] = gauss_integration(A_ii, x[self.index_i], x[self.index_j])
        A_local_[0][1] = gauss_integration(A_ij, x[self.index_i], x[self.index_j])
        A_local_[1][0] = gauss_integration(A_ij, x[self.index_i], x[self.index_j])
        A_local_[1][1] = gauss_integration(A_jj, x[self.index_i], x[self.index_j])
        return A_local_
    
    def M_local(self):
        M_local_ = np.zeros([2,2])
        M_local_[0][0] = gauss_integration(M_ii, x[self.index_i], x[self.index_j])
        M_local_[0][1] = gauss_integration(M_ij, x[self.index_i], x[self.index_j])
        M_local_[1][0] = gauss_integration(M_ij, x[self.index_i], x[self.index_j])
        M_local_[1][1] = gauss_integration(M_jj, x[self.index_i], x[self.index_j])
        return M_local_

    def F_local(self):
        F_local_ = np.zeros(2)
        F_local_[0] = gauss_integration(F_i, x[self.index_i], x[self.index_j])
        F_local_[1] = gauss_integration(F_j, x[self.index_i], x[self.index_j])
        return F_local_
    


# Local element list:
local_list = []
for i in range(dof - 1):
    local_list.append(Local_element(i, i+1))

#2. Stiffness matrix generation:
A = np.zeros([dof, dof])
for local_element in local_list:
    A[local_element.index_i][local_element.index_i] += local_element.A_local()[0][0]
    A[local_element.index_i][local_element.index_j] += local_element.A_local()[0][1]
    A[local_element.index_j][local_element.index_i] += local_element.A_local()[1][0]
    A[local_element.index_j][local_element.index_j] += local_element.A_local()[1][1]

M = np.zeros([dof, dof])
for local_element in local_list:
    M[local_element.index_i][local_element.index_i] += local_element.M_local()[0][0]
    M[local_element.index_i][local_element.index_j] += local_element.M_local()[0][1]
    M[local_element.index_j][local_element.index_i] += local_element.M_local()[1][0]
    M[local_element.index_j][local_element.index_j] += local_element.M_local()[1][1]

A = A + M

#3. Source vector generation:
F = np.zeros(dof)
for local_element in local_list:
    F[local_element.index_i] += local_element.F_local()[0]
    F[local_element.index_j] += local_element.F_local()[1]

#4. Impose boundary conditions:
#impose_Dirichlet(A, F, 0, g1)
#impose_Dirichlet(A, F, dof - 1, g2)
h1 = 0
h2 = 0
impose_Neumann(F, 0, -h1)
impose_Neumann(F, dof - 1, h2)

# 5. Solve the linear system
U = np.linalg.solve(A, F)

# Analytical solution:
#U_an = np.zeros(dof)
#for i in range(dof):
    #U_an[i] = math.sin(x[i])/2
    #U_an[i] = 0.5*(math.sin(x[i]) - math.sinh(x[i]) + math.tanh(math.pi)*math.cosh(x[i]))

# 6. Visualize the solution
plt.plot(x, U , label = "Finite Element Method")
plt.plot(x, U_an,'--r', label = "Analytical Solution")
plt.xlabel('x')
plt.ylabel('y')
plt.title('Finite Element Solver: 1D Laplacian')
plt.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1)
plt.show()