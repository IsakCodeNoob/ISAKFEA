# This document is responsible for loading and parsing input data for the ISAKFEA program.
# It reads input files, extracts necessary parameters, and prepares data structures for analysis.
# The software is developed by Jonas Isaksen, and may be used and modified freely for
# non-commercial purposes, provided that this notice is preserved.

# Packages used in Input_Loader are imported below
import numpy as np                      # Fundamental package for numerical computations
import scipy.sparse as sp               # Sparse matrix package
import re


def EmptyGlobal():
    # Making empty global variables to be used across modules
    G = {}  # Global dictionary to store various parameters and data
    G['NumNodes'] = 0          # Number of nodes in the finite element model
    G['NodeCoords'] = None   # Coordinates of nodes

    G['NumElem'] = 0       # Number of elements in the finite element model
    G['ElemNodes'] = None    # Nodes associated with each element

    G['StiffnessMatrix'] = None # Global stiffness matrix
    G['LoadVector'] = None     # Global load vector
    G['BoundaryConditions'] = []  # Boundary conditions
    G['DisplacementVector'] = None  # Global displacement vector
    G['E'] = 0.0                # Young's modulus [MPa]
    G['t'] = 0.0                # Thickness of elements
    G['nu'] = 0.0               # Poisson's ratio


    return G


# Loading and making input data for the progam

def InitializeGlobal(data,G):
    for line in data:
        Line = line.strip()

        if Line.startswith('NUMOFF,NODE,'): # Number of nodes
            G['NumNodes'] = int(Line.split(',')[2])
            G['NodeCoords'] = np.zeros((G['NumNodes'], 2))  # Initialize node coordinates array
            G['LoadVector'] = np.zeros(2*G['NumNodes'])  # Initialize load vector
            G['DisplacementVector'] = np.zeros(2*G['NumNodes'])  # Initialize displacement vector


        elif Line.startswith('NUMOFF,ELEM,'):   # Number of elements
            G['NumElem'] = int(Line.split(',')[2])
            G['ElemNodes'] = np.zeros((G['NumElem'], 8), dtype=int)  # Initialize element nodes array


# Sub-function for looping through lines of input data
def LineLoops(data, G):
    # Looping through each line to extract relevant information
    for line in data:
        Line = line.strip()

        # Reading node coordinates
        if Line.startswith('(3i9,6e21.13e3)'): # Format specifier for node coordinates
            # Line index for node coordinates starts after the header line
            idx = data.index(line) + 1

            for i in range(G['NumNodes']):

                #if "-" in line:
                #    line.replace("-"," -")

                node_data = data[idx + i]

                # If there is a minus sign, the inp file does not place a space before
                node_data = re.sub(r'(?<!E)-', ' -', node_data)

                node_data = node_data.strip().split(' ')         # Split line into components
                
                
                

                node_data = [float(s) for s in node_data if s != ""] # Convert to float and remove empty strings

                # ANSYS does not write 0.0 in y-coordinate if it is zero, so we handle that case
                if len(node_data) < 5:
                    node_data.insert(4, 0.0)  # Insert 0.0 for Y coordinate if missing

                G['NodeCoords'][i, 0] = node_data[3]  # X coordinate
                G['NodeCoords'][i, 1] = node_data[4]  # Y coordinate

        # Reading element connectivity
        elif Line.startswith('(19i10)'): # Format specifier for element connectivity
            # Line index for element connectivity starts after the header line
            idx = data.index(line) + 1

            for i in range(G['NumElem']): # Looping through each element
                elem_data = data[idx + i].strip().split(' ')         # Split line into components
                elem_data = [int(s) for s in elem_data if s != ""]   # Convert to int and remove empty strings

                G['ElemNodes'][i,:] = elem_data[11:19]  # Store element node indices
        
        # Reading material properties
        elif Line.startswith('MPDATA,UNBL, 1,EX  ,'):   # Young's modulus
            G['E'] = float(Line.split(',')[6])

        elif Line.startswith('MPDATA,UNBL, 1,NUXY,'):   # Poisson's ratio
            G['nu'] = float(Line.split(',')[6])

        elif Line.startswith('(7g16.9)'): # Format specifier for thickness
            # Line index for thickness starts after the header line
            idx = data.index(line) + 1

            thickness_data = data[idx].strip().split(' ')         # Split line into components
            thickness_data = [float(s) for s in thickness_data if s != ""] # Convert to float and remove empty strings

            G['t'] = thickness_data[2]  # Assuming uniform thickness for all elements

        elif Line.startswith('FND'):   # Nodal forces
            
            # removing spaces from each entry
            SplitLine = [entry.strip() for entry in Line.split(',')]
            idx = float(SplitLine[1])  # Node ID

            if SplitLine[2] == 'FX':
                force_value = float(SplitLine[3])
                G['LoadVector'][2*(int(idx)-1)] += force_value  # X direction force
            elif SplitLine[2] == 'FY':
                force_value = float(SplitLine[3])
                G['LoadVector'][2*(int(idx)-1)+1] += force_value  # Y direction force
        
        elif Line.startswith('DND'):    # Displacement boundary conditions
            # removing spaces from each entry
            SplitLine = [entry.strip() for entry in Line.split(',')]
            idx = float(SplitLine[1])  # Node ID

            if SplitLine[2] == 'UX':
                G['BoundaryConditions'].append((2*(int(idx)-1)))  # X direction displacement
            elif SplitLine[2] == 'UY':
                G['BoundaryConditions'].append((2*(int(idx)-1)+1))  # Y direction displacement
            elif SplitLine[2] == 'ALL':
                G['BoundaryConditions'].append((2*(int(idx)-1)))  # X direction displacement
                G['BoundaryConditions'].append((2*(int(idx)-1)+1))  # Y direction displacement

    return G

# Main loading function
def load_input(file_path, G):
    print("Loading input file:", file_path)
    # Loading the data from the input file
    with open(file_path, 'r') as file:   
        data = file.readlines()

    
    # Looping through each line to extract relevant information
    print("Parsing input data...")
    InitializeGlobal(data,G)

    LineLoops(data, G)

    # Making the global stiffness matrix and load vector
    G['StiffnessMatrix'] = sp.lil_matrix((2*G['NumNodes'], 2*G['NumNodes']))  # Initialize global stiffness matrix
    
    return G

