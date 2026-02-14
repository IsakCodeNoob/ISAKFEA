# This module is responsible for computing stress and strain in the ISAKFEA program.
# It calculates strains from displacements and then computes stresses using the constitutive relations.
# The software is developed by Jonas Isaksen, and may be used and modified freely for
# non-commercial purposes, provided that this notice is preserved.

# Packages used in StressStrain are imported below
import numpy as np                      # Fundamental package for numerical computations


# Modules used
import StiffnessMatrix

def ElementDisplacements(Elem, G, elemID):
    

    # Extracting displacement vector for the element
    # Bookkeeping matrix
    LM = np.zeros((2*Elem['NumNodes'],), dtype=int)  # Local to global mapping (element nodeID related to global DOF 2*nodeID-1 and 2*nodeID)
    for i in range(Elem['NumNodes']):
        LM[2*i]   = 2*(Elem['NodeIDs'][i]-1)     # DOF in X direction
        LM[2*i+1] = 2*(Elem['NodeIDs'][i]-1) + 1 # DOF in Y direction
    
    Elem['DisplacementVector'] = G['DisplacementVector'][LM]  # Element displacement vector

    return Elem

def ComputeStressStrain(Elem,C):

    

    # Element strains are computed in the superconvergent points (4 points for PLANE183 elements) [See Cook et al., 2002 Ch. 6.10]
    SuperConvergentPoints = [(-1/np.sqrt(3), -1/np.sqrt(3)), # (xi_1, eta_1)
                             ( 1/np.sqrt(3), -1/np.sqrt(3)), # (xi_2, eta_2)
                             ( 1/np.sqrt(3),  1/np.sqrt(3)), # (xi_3, eta_3)
                             (-1/np.sqrt(3),  1/np.sqrt(3))] # (xi_4, eta_4)
    
    Elem['Strains'] = np.zeros((8,3))  # Initializing strain array for all 8 nodes (epsilon_xx, epsilon_yy, gamma_xy)
    Elem['Stresses'] = np.zeros((8,4))  # Initializing stress array for all 8 nodes (sigma_xx, sigma_yy, tau_xy)
    SuperconStrains = np.zeros((4,3))  # Initializing strain array for superconvergent points
    SuperconStresses = np.zeros((4,3))  # Initializing stress array for superconvergent points

    # Looping over superconvergent points to compute strains and stresses
    for i in range(4):
        J = StiffnessMatrix.ComputeJ(Elem, SuperConvergentPoints[i][0], SuperConvergentPoints[i][1])  # Jacobian matrix at superconvergent point
        J_inv = np.linalg.inv(J)  # Inverse of Jacobian matrix
        B = StiffnessMatrix.ComputeB(Elem, SuperConvergentPoints[i][0], SuperConvergentPoints[i][1], J_inv)  # Strain-displacement matrix at superconvergent point

        SuperconStrains[i,:] = B @ Elem['DisplacementVector']  # Strain at superconvergent point

        SuperconStresses[i,:] = C @ SuperconStrains[i,:]  # Stress at superconvergent point
    

    # Now we need to extrapolate the values to the nodes
    ExtrapolationMatrix = np.sqrt(3)*np.array([
        [-1, -1],
        [ 1, -1],
        [ 1,  1],
        [-1,  1],
        [ 0, -1],
        [ 1,  0],
        [ 0,  1],
        [-1,  0]
    ])

    # Extrapolate strains and stresses to nodes
    for i in range(8):
        N = np.zeros((4)) # Shape functions at superconvergent points
        N = StiffnessMatrix.ComputeNExtrap(N, ExtrapolationMatrix[i,0], ExtrapolationMatrix[i,1])

        # Looping over stress and strain components to extrapolate
        for j in range(3):
            Elem['Strains'][i,j] = N @ SuperconStrains[:,j]
            Elem['Stresses'][i,j] = N @ SuperconStresses[:,j]

        # von Mises equivalent stress
        Elem['Stresses'][i,3] = np.sqrt(
                                        Elem['Stresses'][i,0]**2 -                              # sigma_xx^2
                                        Elem['Stresses'][i,0]*Elem['Stresses'][i,1] +           # - sigma_xx * sigma_yy
                                        Elem['Stresses'][i,1]**2 + 3*Elem['Stresses'][i,2]**2   # + sigma_yy^2 + 3*tau_xy^2
                                        )
        
    return Elem

