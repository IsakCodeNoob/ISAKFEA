# This document is responsible for solving the system of equations for the ISAKFEA program.
# It solves the system of equations K * D = P, where K is the global stiffness matrix,
# D is the displacement vector, and P is the force vector.
# The software is developed by Jonas Isaksen, and may be used and modified freely for
# non-commercial purposes, provided that this notice is preserved.


# Packages used in Sol are imported below
import numpy as np                      # Fundamental package for numerical computations
import scipy.sparse.linalg as spla      # Sparse linear algebra package


def SystemSetup(G):
    # Setting up the system of equations by applying boundary conditions

    total_DOF = G['NumNodes'] * 2  # Total degrees of freedom (2 DOF per node in 2D problems)

    c = np.array(G['BoundaryConditions'])  # Constrained DOF from boundary conditions
    u = np.setdiff1d(np.arange(total_DOF), c)  # Unconstrained DOF

    # Extracting sub-matrices and vectors for unconstrained DOF
    # Assuming that the applied BCs are zero displacements
    Kuu = G['StiffnessMatrix'][u,:][:,u]  # Stiffness matrix for unconstrained DOF
    Pu = G['LoadVector'][u]                # Load vector for unconstrained DOF
    G['Kuu'] = Kuu
    G['Pu'] = Pu

    return G


def SolveSystem(G):
    # Solving the system of equations for unconstrained DOF

    # Using sparse linear solver for efficiency
    Du = spla.spsolve(G['Kuu'], G['Pu'])  # Displacement vector for unconstrained DOF

    G['DisplacementVectorFree'] = Du  # Storing displacement vector in global dictionary

    return G


def AssembleDisplacementVector(G):
    # Assembling the full displacement vector including constrained DOF

    total_DOF = G['NumNodes'] * 2  # Total degrees of freedom (2 DOF per node in 2D problems)

    c = np.array(G['BoundaryConditions'])  # Constrained DOF from boundary conditions
    u = np.setdiff1d(np.arange(total_DOF), c)  # Unconstrained DOF

    # Assigning computed displacements to unconstrained DOF
    G['DisplacementVector'][u] = G['DisplacementVectorFree']

    return G

