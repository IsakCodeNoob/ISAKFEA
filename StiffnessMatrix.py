# This document is responsible for constructing the global stiffness matrix for the ISAKFEA program.
# It assembles the stiffness contributions from individual elements into a global matrix.
# The only considered element type is PLANE183 elements (8-nnode isoparametric quadrilateral elements).
# The software is developed by Jonas Isaksen, and may be used and modified freely for
# non-commercial purposes, provided that this notice is preserved.

# Packages used in StiffnessMatrix are imported below
import numpy as np                      # Fundamental package for numerical computations



# The empty local element dictionary (Assumption of PLANE183 elements)
Elem = {}  # Local dictionary to store element-specific parameters and data
def Empty_Local(G,elem_id):

    Elem['NumNodes'] = 8        # Number of nodes per element ( PLANE183 element )
    Elem['Ke'] = np.zeros((2*Elem['NumNodes'], 2*Elem['NumNodes']))           # Element stiffness matrix
    Elem['NodeIDs'] = None      # Node IDs associated with the element (used for bookkeeping)
    Elem['NodeIDs'] = G['ElemNodes'][elem_id,:]  # Extracting node IDs for the element
    Elem['ElemID'] = elem_id    # Element ID

    Elem['Coords'] = None       # Coordinates of the element nodes
    Elem['Coords'] = G['NodeCoords'][Elem['NodeIDs']-1,:]  # Extracting node coordinates for the element (-1 for zero indexing done in python)

    Elem['E'] = G['E']          # Young's modulus [MPa]
    Elem['t'] = G['t']          # Thickness of elements
    Elem['nu'] = G['nu']        # Poisson's ratio
    

    return Elem


def ComputeJ(Elem,xi,eta):
    # Computing the Jacobian matrix for the element at given natural coordinates (xi, eta)

    # Shape function derivatives with respect to natural coordinates
    dN_nat = np.zeros((2,Elem['NumNodes']))

    # Derivatives w.r.t. xi
    dN_nat[0,0] = -0.25*(eta+2*xi)*(eta-1)
    dN_nat[0,1] = 0.25*(eta-2*xi)*(eta-1)
    dN_nat[0,2] = 0.25*(eta+2*xi)*(eta+1)
    dN_nat[0,3] = -0.25*(eta-2*xi)*(eta+1)
    dN_nat[0,4] =  xi*(eta-1)
    dN_nat[0,5] = -0.5*eta**2+0.5
    dN_nat[0,6] = -xi*(eta+1)
    dN_nat[0,7] = 0.5*eta**2-0.5

    # Derivatives w.r.t. eta
    dN_nat[1,0] = -0.25*(xi+2*eta)*(xi-1)
    dN_nat[1,1] = -0.25*(xi-2*eta)*(xi+1)
    dN_nat[1,2] = 0.25*(xi+2*eta)*(xi+1)
    dN_nat[1,3] = 0.25*(xi-2*eta)*(xi-1)
    dN_nat[1,4] = 0.5*xi**2-0.5
    dN_nat[1,5] = -eta*(xi+1)
    dN_nat[1,6] = -0.5*xi**2+0.5
    dN_nat[1,7] = eta*(xi-1)

    # Initializing Jacobian matrix
    J = np.zeros((2,2))

    # Computing Jacobian matrix components
    Coords = Elem['Coords']
    
    J = dN_nat @ Coords

    return J

def ComputeB(Elem,xi,eta,J_inv):
    # Computing the strain-displacement matrix B for the element at given natural coordinates (xi, eta)

    # Shape function derivatives with respect to natural coordinates
    dN_nat = np.zeros((2,Elem['NumNodes']))

    # Derivatives w.r.t. xi
    dN_nat[0,0] = -0.25*(eta+2*xi)*(eta-1)
    dN_nat[0,1] = 0.25*(eta-2*xi)*(eta-1)
    dN_nat[0,2] = 0.25*(eta+2*xi)*(eta+1)
    dN_nat[0,3] = -0.25*(eta-2*xi)*(eta+1)
    dN_nat[0,4] =  xi*(eta-1)
    dN_nat[0,5] = -0.5*eta**2+0.5
    dN_nat[0,6] = -xi*(eta+1)
    dN_nat[0,7] = 0.5*eta**2-0.5

    # Derivatives w.r.t. eta
    dN_nat[1,0] = -0.25*(xi+2*eta)*(xi-1)
    dN_nat[1,1] = -0.25*(xi-2*eta)*(xi+1)
    dN_nat[1,2] = 0.25*(xi+2*eta)*(xi+1)
    dN_nat[1,3] = 0.25*(xi-2*eta)*(xi-1)
    dN_nat[1,4] = 0.5*xi**2-0.5
    dN_nat[1,5] = -eta*(xi+1)
    dN_nat[1,6] = -0.5*xi**2+0.5
    dN_nat[1,7] = eta*(xi-1)

    # Shape function derivatives with respect to global coordinates
    dN_glob = np.zeros((2,Elem['NumNodes']))

    for ii in range(Elem['NumNodes']):
        for jj in range(2):
            dN_glob[jj,ii] = J_inv[jj,0]*dN_nat[0,ii] + J_inv[jj,1]*dN_nat[1,ii]

    # Constructing the B matrix
    B = np.zeros((3, 2*Elem['NumNodes']))
    for k in range(Elem['NumNodes']):
        B[0, 2*k]     = dN_glob[0,k]
        B[1, 2*k + 1] = dN_glob[1,k]
        B[2, 2*k]     = dN_glob[1,k]
        B[2, 2*k + 1] = dN_glob[0,k]

    return B


def ComputeNExtrap(N,r,s):
    
    # Shape functions for 4-node superconvergent points
    N[0] = 0.25*(1-r)*(1-s)
    N[1] = 0.25*(1+r)*(1-s)
    N[2] = 0.25*(1+r)*(1+s)
    N[3] = 0.25*(1-r)*(1+s)

    return N


def LocalStiffness(Elem):

    # Material properties
    E = Elem['E']
    nu = Elem['nu']
    t = Elem['t']

    # Constitutive matrix for plane stress
    C = np.array([
        [E/(1-nu**2),  nu*E/(1-nu**2), 0],
        [nu*E/(1-nu**2), E/(1-nu**2), 0],
        [0, 0, E/(2*(1+nu))]
    ])

    # Integration points (Gauss quadrature, reduced integration)
    xi =  1/np.sqrt(3)*np.array([-1, 1])
    eta = 1/np.sqrt(3)*np.array([-1, 1])

    # Gauss integration to compute element stiffness matrix
    for i in range(len(xi)):
        for j in range(len(eta)):
            J = ComputeJ(Elem, xi[i], eta[j])
            J_inv = np.linalg.inv(J)
            B = ComputeB(Elem, xi[i], eta[j], J_inv)

            Elem['Ke'] += t * np.linalg.det(J) * (B.T @ C @ B)

    return Elem['Ke']