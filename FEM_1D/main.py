# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE FEM 1D
Title: main.py
"""
# Basic packages:
from config_packages import np, math, plt, cm, logging, data, module_name

# Configure the logging format and level (adjust as needed)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")

# Custom packages:
from src import mesh_generation, boundary_conditions, post_processing

def main():
    #0. Read data:
    logging.info(f"(0/6) Reading data from {module_name}.py ...")

    #1. Mesh generation:
    logging.info("(1/6) Generating mesh ...")
    mesh = mesh_generation.Mesh([data.x1, data.x2], data.h)

    #2. Assemble of global matrices and right hand side:
    logging.info("(2/6) Assembling global matrices and right hand side ...")
    A = np.zeros([mesh.ndof, mesh.ndof])
    element = 0 
    for local_element in mesh.elements:
        for m in range(2):
            for n in range(2):
                A[mesh.elements_nodes[element][m]][mesh.elements_nodes[element][n]] += local_element.A_local()[m][n]
        element = element + 1 

    F = np.zeros(mesh.ndof)
    element = 0
    for local_element in mesh.elements:
        for m in range(2):
            F[mesh.elements_nodes[element][m]] += local_element.F_local()[m]
        element = element + 1

    # 3. Impose boundary conditions:
    logging.info("(3/6) Imposing boundary conditions ...")
    A, F, g = boundary_conditions.impose_boundary_conditions(A, F, mesh)

    # 4. Solve the algebraic problem:
    logging.info("(4/6) Solving the algebric problem ...")
    U = np.linalg.solve(A, F)
    # Lifting operation:
    U = U + g

    # 5. Plotting the solution:
    logging.info("(5/6) Plotting the solution ...")
    # Plot the solution
    if data.plot_solution == 'y':
        fig = plt.figure()
        plt.plot(mesh.coord, U , label = "Finite Element Method")
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Finite Element Solver: 1D Laplacian')
        plt.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1)
        plt.draw()

    # 6. Computing the error:
    #logging.info("(6/6) Computing errors ...")
    #err = post_processing.compute_errors(U, mesh)

    # 7. Finish the program:
    plt.show()
    exit()

if __name__ == '__main__':
    main()