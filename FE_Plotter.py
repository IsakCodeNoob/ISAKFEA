# This is plotting document for the ISAKFEA program
# The software is developed by Jonas Isaksen, and may be used and modified freely for
# non-commercial purposes, provided that this notice is preserved.

# Packages used in ISAKFEA are imported below
import numpy as np                      # Fundamental package for numerical computations
import matplotlib.pyplot as plt         # Package for plotting results
import matplotlib.colors as mcolors     # Package for handling colors in plots
from tqdm import tqdm                   # Progressbar

# Importing core modules of ISAKFEA
import Input_Loader

# Importing input data using the Input_Loader module
G = Input_Loader.EmptyGlobal()  # Initialize global variables

file_name = input('Enter the name of the input file: ')
file_path = f"IO/{file_name}.inp"  # Constructing file path
G = Input_Loader.load_input(file_path, G)  # Module for loading and parsing input data
print("Input data loaded successfully.")

# Importing results from output file
results_file_path = f"IO/{file_name}/{file_name}_"  # Constructing file path for results (assuming results are stored in IO/{file_name}/{file_name}_*.txt)

# Determine the scale factor for plotting displacements
DispScale = input("Enter the scale factor for plotting displacements (e.g., 0 for no scaling, 10 for 10x scaling): ")
# Loading displacement files
Nux = np.loadtxt(results_file_path + "Nux.txt",delimiter=',',skiprows=1)[:,1]  # Load displacement in x direction
Nuy = np.loadtxt(results_file_path + "Nuy.txt",delimiter=',',skiprows=1)[:,1]  # Load displacement in y direction
# Scaling displacements for plotting
Nux = Nux*float(DispScale)  # Scale displacement in x direction
Nuy = Nuy*float(DispScale)  # Scale displacement in y direction
# Original coordinates of nodes + scaled displacements
X = G['NodeCoords'][:,0] + Nux  # X coordinates of nodes
Y = G['NodeCoords'][:,1] + Nuy  # Y coordinates of nodes


# User input for which results to plot
print("Select the type of plot to display:")
print("1: Displacement in x direction")
print("2: Displacement in y direction")
print("3: Stress in xx direction")
print("4: Stress in yy direction")
print("5: Stress in xy direction")
print("6: Von Mises stress")
print("7: Strain in xx direction")
print("8: Strain in yy direction")
print("9: Shear strain in xy direction")
plot_type = float(input("Enter the type of plot (displacement/stress/strain): "))
# Loading the appropriate results based on user input
if plot_type == 1:
    # Load displacement results for x
    PlotVar = np.loadtxt(results_file_path + "Nux.txt",delimiter=',',skiprows=1)[:,1]  # Load displacement in x direction
    title = 'FEA x-displacement'
    ColorbarLabel = 'Displacement in x direction [mm]'
elif plot_type == 2:
    # Load displacement results for y
    PlotVar = np.loadtxt(results_file_path + "Nuy.txt",delimiter=',',skiprows=1)[:,1]  # Load displacement in y direction
    title = 'FEA y-displacement'
    ColorbarLabel = 'Displacement in y direction [mm]'
elif plot_type == 3:
    # Looad stress results for xx
    PlotVar = np.loadtxt(results_file_path + "Nsx.txt",delimiter=',',skiprows=1)[:,1]  # Load stress in xx direction
    title = 'FEA x-stress'
    ColorbarLabel = 'Stress in xx direction [MPa]'
elif plot_type == 4:
    # Looad stress results for yy
    PlotVar = np.loadtxt(results_file_path + "Nsy.txt",delimiter=',',skiprows=1)[:,1]  # Load stress in yy direction
    title = 'FEA y-stress'
    ColorbarLabel = 'Stress in yy direction [MPa]'
elif plot_type == 5:
    # Looad stress results for xy
    PlotVar = np.loadtxt(results_file_path + "Ntxy.txt",delimiter=',',skiprows=1)[:,1]  # Load stress in xy direction
    title = 'FEA xy-stress'    
    ColorbarLabel = 'Shear stress in xy direction [MPa]'
elif plot_type == 6:
    # Looad stress results for NsRef
    PlotVar = np.loadtxt(results_file_path + "NsRef.txt",delimiter=',',skiprows=1)[:,1]  # Load von mises stress
    title = 'FEA von Mises stress'
    ColorbarLabel = 'Von Mises stress [MPa]'
elif plot_type == 7:
    # Looad strain results for xx
    PlotVar = np.loadtxt(results_file_path + "Nepsx.txt",delimiter=',',skiprows=1)[:,1]*1e6  # Load strain in xx direction
    title = 'FEA x-strain'
    ColorbarLabel = 'Strain in xx direction [mum/m]'
elif plot_type == 8:
    # Looad strain results for yy
    PlotVar = np.loadtxt(results_file_path + "Nepsy.txt",delimiter=',',skiprows=1)[:,1]*1e6  # Load strain in yy direction
    title = 'FEA y-strain'
    ColorbarLabel = 'Strain in yy direction [mum/m]'
elif plot_type == 9:
    # Looad strain results for xy
    PlotVar = np.loadtxt(results_file_path + "Ngamxy.txt",delimiter=',',skiprows=1)[:,1]*1e6  # Load shear strain in xy direction
    title = 'FEA xy-strain'    
    ColorbarLabel = 'Shear strain in xy direction [mum/m]'
else:
    print("Invalid input. Please enter a number between 1 and 9.")
    exit()




# Plotting results

# Defining levels for the plot
CustomLevels = input("Press 1 for custom levels or press Enter to use automatic levels: ")
if CustomLevels == "1":
    print("Current range of values: ", min(PlotVar), " to ", max(PlotVar))
    min_level = float(input("Enter second lowest level: "))
    max_level = float(input("Enter second highest level: "))
    num_levels = int(input("Enter number of levels: "))
    ContourLevels = np.linspace(min_level, max_level, num_levels - 2)  # Generate levels between min and max
    ContourLevels = np.concatenate(([min(PlotVar)], ContourLevels, [max(PlotVar)]))  # Add min and max levels to the array
else:
    ContourLevels = np.linspace(min(PlotVar), max(PlotVar), 20)

plt.figure()  # Set figure size

# The following two lines create a colormap and normalization based on the contour levels, ensuring that the colors in the plot correspond to the defined levels.
cmap = plt.get_cmap('jet', len(ContourLevels)-1)  # Get the 'jet' colormap with the number of levels
norm = mcolors.BoundaryNorm(ContourLevels, cmap.N)  # Create a normalization object based on the contour levels

# Looping over each element in the model
for ii in tqdm(range(G['NumElem'])):
    
    nodeidx = G['ElemNodes'][ii,0:4] - 1 # so far only plotting corner nodes (-1 due to zero indexing)

    #ElemX , ElemY = np.meshgrid(X[nodeidx], Y[nodeidx])
    ElemX = X[nodeidx]
    ElemY = Y[nodeidx]
    ElemVar = PlotVar[nodeidx]

    plt.tricontourf(ElemX,ElemY,ElemVar,cmap=cmap,levels=ContourLevels,norm=norm)  # Create filled contour plot for the element


plt.colorbar(label=ColorbarLabel)  # Add colorbar with label
plt.title(title)  # Set title
plt.xlabel('X Coordinate')  # Set x-axis label
plt.ylabel('Y Coordinate')  # Set y-axis label
plt.axis('equal')  # Set equal aspect ratio for x and y axes
plt.ylim([min(Y)-0.1*max(Y), max(Y)+0.1*max(Y)])  # Set y-axis limits with some padding
plt.xlim([min(X)-0.1*max(X), max(X)+0.1*max(X)])  # Set x-axis limits with some padding
plt.show()  # Display the plot


