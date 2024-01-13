# Author: Cássio Murakami
# Finite Element Solver: 1D Laplacian
import numpy as np 
import math
import matplotlib.pyplot as plt

def f(x):
    return math.sin(x)

def impose_Dirichlet(A, F, index, g):
    A_ = A
    F_ = F
    for i in range(len(A)):
        A_[index][i] = 0
    A_[index][index] = 1
    F_[index] = g
    return (A_, F_)

def gauss_integration(basis_function, a, b):
    # Gaussian quadrature parameters:
    zeta_list = [-math.sqrt(3/5), 0 , + math.sqrt(3/5)]
    weight_list = [5/9, 8/9 , 5/9]

    # Gaussian integration loop:
    result = 0
    for q in range(len(zeta_list)):
        x_ = zeta_list[q]*(b - a)/2 + (a + b)/2
        if basis_function == 1:
            phi = (x_ - a)/(b - a)
        elif basis_function == 2: 
            phi = (b - x_)/(b - a)

        result = result + 0.5*(b-a)*weight_list[q]*phi*f(x_)
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
        # Definition of the local stiffness matrix:
        A_local_ = np.zeros([2, 2])
        A_local_[0][0] = 1/(x[self.index_j] - x[self.index_i])
        A_local_[0][1] = -1/(x[self.index_j] - x[self.index_i])
        A_local_[1][0] = -1/(x[self.index_j] - x[self.index_i])
        A_local_[1][1] = 1/(x[self.index_j] - x[self.index_i])
        return A_local_
    
    def F_local(self):
        # Definition of the local source vector:
        F_local_ = np.zeros(2)
        F_local_[0] = gauss_integration(1, x[self.index_i], x[self.index_j])
        F_local_[1] = gauss_integration(2, x[self.index_i], x[self.index_j])
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

#3. Source vector generation:
F = np.zeros(dof)
for local_element in local_list:
    F[local_element.index_i] += local_element.F_local()[0]
    F[local_element.index_j] += local_element.F_local()[1]

#4. Impose boundary conditions:
impose_Dirichlet(A, F, 0, g1)
impose_Dirichlet(A, F, dof - 1, g2)

# 5. Solve the linear system
U = np.linalg.solve(A, F)

# Analytical solution:
U_an = np.zeros(dof)
for i in range(dof):
    U_an[i] = math.sin(x[i])

    
# 6. Visualize the solution
plt.plot(x, U , label = "Finite Element Method")
plt.plot(x, U_an,'--r', label = "Analytical Solution")
plt.xlabel('x')
plt.ylabel('y')
plt.title('Finite Element Solver: 1D Laplacian')
plt.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1)
plt.show()