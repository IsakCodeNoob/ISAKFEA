# This is the main document of the ISAKFEA program.
# ISAKFEA is a finite element analysis software for structural engineering applications.
# For more information, read the README file included in the distribution package.
# The software is developed by Jonas Isaksen, and may be used and modified freely for
# non-commercial purposes, provided that this notice is preserved.

# Packages used in ISAKFEA are imported below
import numpy as np                      # Fundamental package for numerical computations
import scipy.sparse as sp               # Sparse matrix package
import scipy.sparse.linalg as spla      # Sparse linear algebra package
import timeit                           # Package for measuring execution time
from tqdm import tqdm                   # Progressbar

# Importing core modules of ISAKFEA
import Input_Loader
import StiffnessMatrix
import Sol
import StressStrain
import OutputFiles



# Importing input data using the Input_Loader module
G = Input_Loader.EmptyGlobal()  # Initialize global variables

file_name = input('Enter the name of the input file: ')
file_path = f"IO/{file_name}.inp"  # Constructing file path
G = Input_Loader.load_input(file_path, G)  # Module for loading and parsing input data
print("Input data loaded successfully.")

# Constructing stiffness matrix using the StiffnessMatrix module
# Global stiffness matrix is assembled by K = LeT * Ke * Le for each element
print("Constructing global stiffness matrix...")
Starttime = timeit.default_timer()  # Start timer
for elem_id in tqdm(range(G['NumElem'])):
    Elem = StiffnessMatrix.Empty_Local(G,elem_id)  # Initialize local element dictionary
    Elem['Ke'] = StiffnessMatrix.LocalStiffness(Elem)  # Compute local stiffness matrix

    # Assembling the bookkeeping matrix
    LM = np.zeros((2*Elem['NumNodes'],), dtype=int)  # Local to global mapping (element nodeID related to global DOF 2*nodeID-1 and 2*nodeID)
    for i in range(Elem['NumNodes']):
        LM[2*i]   = 2*(Elem['NodeIDs'][i]-1)     # DOF in X direction
        LM[2*i+1] = 2*(Elem['NodeIDs'][i]-1) + 1 # DOF in Y direction

    L = sp.lil_matrix((2*Elem['NumNodes'], G['NumNodes']*2))  # Bookkeeping matrix (16 x total DOF)
    # Filling the bookkeeping matrix (ones at appropriate locations)
    for i in range(2*Elem['NumNodes']):
        L[i,LM[i]] = 1.0

    # Assembling global stiffness matrix
    G['StiffnessMatrix'] += L.T @ Elem['Ke'] @ L  # K = LeT * Ke * Le (Cook eq 4.8-18 to 4.8-20)

# Inefficient form of sparce matrix. Changing to CSR format after assembly for efficiency (better for solving)
G['StiffnessMatrix'] = G['StiffnessMatrix'].tocsr()  # Convert to CSR format for efficiency


endtime = timeit.default_timer()    # End timer
print("Global stiffness matrix constructed successfully.")
print(f"Global stiffness matrix constructed in: {endtime - Starttime} seconds.")


print("Setting up system of equations...")
G = Sol.SystemSetup(G)
print("System of equations set up successfully.")

# Storing global stiffnes matrix is not necessary anymore, freeing up memory
G['StiffnessMatrix'] = None


print("Solving system of equations...")
Starttime = timeit.default_timer()  # Start timer
G = Sol.SolveSystem(G)
endtime = timeit.default_timer()    # End timer
print("System of equations solved successfully.")
print(f"System of equations solved in: {endtime - Starttime} seconds.")

print("Assembling global displacement vector...")
G = Sol.AssembleDisplacementVector(G)
G['DisplacementVectorFree'] = None  # Freeing up memory
print("Global displacement vector assembled successfully.")

print("Computing stresses and strains...")
Starttime = timeit.default_timer()  # Start timer
# Looping through each element to compute stresses and strains

G['Epsxx'] = np.zeros((G['NumNodes'],))  # Initializing global epsilon_xx array
G['Epsyy'] = np.zeros((G['NumNodes'],))  # Initializing global epsilon_yy array
G['Gammaxy'] = np.zeros((G['NumNodes'],)) # Initializing global gamma_xy array
G['Sigmaxx'] = np.zeros((G['NumNodes'],))  # Initializing global sigma_xx array
G['Sigmayy'] = np.zeros((G['NumNodes'],))  # Initializing global sigma_yy array
G['Tauxy'] = np.zeros((G['NumNodes'],))    # Initializing global tau_xy array
G['SigVM'] = np.zeros((G['NumNodes'],))     # Initializing global von Mises stress array

# As multiple element share nodes, we will take the average values. Thus, we will also keep track of the number of contributions per node.
NodeContributions = np.zeros((G['NumNodes'],), dtype=int)  # Initializing node contributions array

for elem_id in range(G['NumElem']):

    Elem = StiffnessMatrix.Empty_Local(G, elem_id)  # Initialize local element dictionary
    Elem = StressStrain.ElementDisplacements(Elem, G, elem_id)  # Extract element displacements

    # Constitutive matrix for plane stress
    C = np.array([
        [Elem['E']/(1-Elem['nu']**2),  Elem['nu']*Elem['E']/(1-Elem['nu']**2), 0],
        [Elem['nu']*Elem['E']/(1-Elem['nu']**2), Elem['E']/(1-Elem['nu']**2), 0],
        [0, 0, Elem['E']/(2*(1+Elem['nu']))]
    ])

    Elem = StressStrain.ComputeStressStrain(Elem, C)  # Compute stresses and strains for the element

    # Assemling into global arrays
    for i in range(Elem['NumNodes']):
        node_id = Elem['NodeIDs'][i] - 1  # Global node ID (0-indexed)
        NodeContributions[node_id] += 1   # Incrementing contribution count for the node
        G['Epsxx'][node_id]  += Elem['Strains'][i,0]
        G['Epsyy'][node_id]  += Elem['Strains'][i,1]
        G['Gammaxy'][node_id] += Elem['Strains'][i,2]

        G['Sigmaxx'][node_id] += Elem['Stresses'][i,0]
        G['Sigmayy'][node_id] += Elem['Stresses'][i,1]
        G['Tauxy'][node_id]   += Elem['Stresses'][i,2]
        G['SigVM'][node_id]    += Elem['Stresses'][i,3]
    
# Averaging the contributions for each node
for node_id in range(G['NumNodes']):
    if NodeContributions[node_id] > 0:
        G['Epsxx'][node_id]  /= NodeContributions[node_id]
        G['Epsyy'][node_id]  /= NodeContributions[node_id]
        G['Gammaxy'][node_id] /= NodeContributions[node_id]

        G['Sigmaxx'][node_id] /= NodeContributions[node_id]
        G['Sigmayy'][node_id] /= NodeContributions[node_id]
        G['Tauxy'][node_id]   /= NodeContributions[node_id]
        G['SigVM'][node_id]    /= NodeContributions[node_id]

C = None # Freeing up memory

endtime = timeit.default_timer()    # End timer
print("Stresses and strains computed successfully.")
print(f"Stresses and strains computed in: {endtime - Starttime} seconds.")

print("Writing results to output files...")
OutputFiles.WriteResults(G, file_name)
print("Results written successfully.")

