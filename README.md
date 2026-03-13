# 2D FEM Transient Heat Transfer Solver

A custom **C++** implementation of the Finite Element Method (FEM) for solving 2D transient heat transfer problems. Built entirely from scratch without external mathematical libraries.

## Key Features
* **Custom Mesh Parser:** Reads and processes nodes, 4-node quadrilateral elements, and boundary conditions directly from structured text files.
* **Numerical Integration:** Uses Gauss quadrature (configurable integration points) for calculating Jacobians, shape functions, and element matrices.
* **Global Matrix Assembly:** Efficiently builds global stiffness (H) and heat capacity (C) matrices, applying Neumann/Dirichlet boundary conditions.
* **Transient Analysis:** Simulates temperature distribution over time using the Euler backward method, tracking Min/Max node temperatures per time step.
* **Custom Math Engine:** Includes a standalone linear system solver using Gaussian elimination with partial pivoting.

## Quick Start
```bash

g++ -std=c++20 main.cpp src/femCalculator.cpp src/mesh.cpp -o fem_solver

./fem_solver
```
