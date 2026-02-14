The ISAKFEA program is a finite element program written by Jonas Isaksen. It may be used and modified freely, 
for non-commercial purposes, provided that this notice is preserved.

The program is written in Python, and zero attempts have been made to optimize the construction of the global stiffness matrix. Thus, the program is very inefficient!

Units for the program is mm, MPa, m/m (for strain)


Main document: ISAKFEA.py
- This document is used as the main document that calls the rest of the modules. Some operations are done in this document, to avoid exessive functions

Modules and functions:

- Input_Loader.py
-- EmptyGlobal: Initializes the empty global variables for the analysis
-- InitializeGlobal: Fills in the number of nodes and numers of elements in the input file. Thus initializing the DOF dependent variables
-- LineLoops: Used to search each line for relevant information. See bottom of README for requirements of the input file such that the function finds the information
-- load_input: Main function in this module, which calls the previous functions

- StiffnessMatrix.py
-- Empty_Local: Initializes local variables used for each element when calculating the elementwise stuff (Stiffness and stress)
-- ComputeJ: Computes Jacobian matrix in defined spots. Thus it is used in loops to calculate Jacobian in e.g. Gauss-Points
-- ComputeB: Computes strain-displacement matrix in defined spots. Thus it is used in loops to calculate B-matrix in e.g. Gauss-Points
-- ComputeNExtrap: Computes shape-functions of a 4-node element in defined points. This is used for stress-extrapolation from superconvergent points
-- LocalStiffness: Main function of this module. Calls ComputeJ and ComputeB in 2x2 gauss points to calculate the element stiffness matrix through reduced integration.

- Sol.py
-- SystemSetup: Used to extract submatrices from the total linear system. We want to solve the system for displacements of the unconstrained nodes. Constrained DOF are assumed to be zero
-- SolveSystem: Solves the system of equations using spla.spsolve. This is only for the unconstrained nodes. 
-- AssembleDisplacementVector: Assembles the total displacement vector by assigning the solved displacements to the unconstrained DOFs

- StressStrain.py
-- ElementDisplacements: From global displacement vector, the element displacements are extracted for a given element ID
-- ComputeStressStrain: Calculates stresses and strains from element displacements. This is done from 2x2 superconvergent points and extrapolation

- OutputFiles.py
-- WriteResults: Creates a folder with the project name in which the solutions are written for each node in the model.


Illustrations (FE_Plotter.py):
- The plotting program for the solutions of ISAKFEA.py
- Currently, only corner nodes for each element are integrated in the plotter. Thus, midside nodes are not plotted


Requirements for input file:
- Filetype: '.inp'
-- Do this through ANSYS mechanical
- Format specifier for node coordinates: '(3i9,6e21.13e3)'
- Format specifier for element connectivity: '(19i10)'
- Format specifier for Young's modulus: 'MPDATA,UNBL, 1,EX  ,'
- Format specifier for Poisson's ratio: 'MPDATA,UNBL, 1,NUXY,'
- Format specifier for model thickness: '(7g16.9)'
- Format specifier for applied forces: 'FND'
- Format specifier for displacement BCs: 'DND'
- For examples and better illustrations, see '3PB.inp', as the program is based on this file. It is a three point bending model.
- Input files must be located in the folder 'IO\'
-- Necessary as the output folder and files are created here. Furthermore, the FE_Plotter.py assumes result files are located here. 



