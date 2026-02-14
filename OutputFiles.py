# This module is responsible for exporting all the calculated results from the ISAKFEA program.
# It writes the nodal displacements, stresses, and strains to output files for further analysis
# The software is developed by Jonas Isaksen, and may be used and modified freely for
# non-commercial purposes, provided that this notice is preserved.

# Packages used in OutputFiles are imported below
import os                               # Operating system interface package

def WriteResults(G, filename):
    # Ensure output directory exists
    output_dir = f"IO/{filename}/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    # Writing nodal displacements to output files
    disp_file = f"IO/{filename}/{filename}_Nux.txt"
    with open(disp_file, 'w') as f:
        f.write("NodeID, Ux\n")
        for node_id in range(G['NumNodes']):
            Ux = G['DisplacementVector'][2*node_id]
            f.write(f"{node_id+1}, {Ux:.6e}\n")

    
    disp_file = f"IO/{filename}/{filename}_Nuy.txt"
    with open(disp_file, 'w') as f:
        f.write("NodeID, Uy\n")
        for node_id in range(G['NumNodes']):
            Uy = G['DisplacementVector'][2*node_id + 1]
            f.write(f"{node_id+1}, {Uy:.6e}\n")
    

    # Writing nodal strains to output files
    strain_file = f"IO/{filename}/{filename}_Nepsx.txt"
    with open(strain_file, 'w') as f:
        f.write("NodeID, Epsx\n")
        for node_id in range(G['NumNodes']):
            Epsx = G['Epsxx'][node_id] 
            f.write(f"{node_id+1}, {Epsx:.6e}\n")
    
    
    strain_file = f"IO/{filename}/{filename}_Nepsy.txt"
    with open(strain_file, 'w') as f:
        f.write("NodeID, Epsy\n")
        for node_id in range(G['NumNodes']):
            Epsy = G['Epsyy'][node_id] 
            f.write(f"{node_id+1}, {Epsy:.6e}\n")
    
    
    strain_file = f"IO/{filename}/{filename}_Ngamxy.txt"
    with open(strain_file, 'w') as f:
        f.write("NodeID, Gamxy\n")
        for node_id in range(G['NumNodes']):
            Gamxy = G['Gammaxy'][node_id] 
            f.write(f"{node_id+1}, {Gamxy:.6e}\n")

    # Writing nodal stresses to output files
    stress_file = f"IO/{filename}/{filename}_Nsx.txt"
    with open(stress_file, 'w') as f:
        f.write("NodeID, Sx\n")
        for node_id in range(G['NumNodes']):
            Sx = G['Sigmaxx'][node_id] 
            f.write(f"{node_id+1}, {Sx:.6e}\n")
    
    
    stress_file = f"IO/{filename}/{filename}_Nsy.txt"
    with open(stress_file, 'w') as f:
        f.write("NodeID, Sy\n")
        for node_id in range(G['NumNodes']):
            Sy = G['Sigmayy'][node_id] 
            f.write(f"{node_id+1}, {Sy:.6e}\n")
    
    
    stress_file = f"IO/{filename}/{filename}_Ntxy.txt"
    with open(stress_file, 'w') as f:
        f.write("NodeID, Tauxy\n")
        for node_id in range(G['NumNodes']):
            Tauxy = G['Tauxy'][node_id] 
            f.write(f"{node_id+1}, {Tauxy:.6e}\n")
    
    
    stress_file = f"IO/{filename}/{filename}_NsRef.txt"
    with open(stress_file, 'w') as f:
        f.write("NodeID, SigVM\n")
        for node_id in range(G['NumNodes']):
            SigVM = G['SigVM'][node_id] 
            f.write(f"{node_id+1}, {SigVM:.6e}\n")
    
    print(f"Results written to IO/{filename}/ successfully.")



