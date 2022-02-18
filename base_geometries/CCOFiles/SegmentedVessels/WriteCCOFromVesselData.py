import numpy as np
import glob
import csv 

FOV    = 0.3                    # Size of the field of view. Shifts the FAZ to (x,y,z)=(0,0,0)
Radius = 0.004                  # Fixed radius for all vessels
                                # TODO: add radius for each vessel
NVessels = 0                    # Number of vessels
NSegments = 0                   # Total number of segments in the tree
Vessels  = []                   # Container for vessels [(start_node, end_node)]

nodes = []                      # Container for nodes [(x0,y0,z0), (x1,y1,z1)...]
connectivity = []               # Container for connectivity 

# Read data from files
for filename in glob.iglob('./Vessel*'):
    NVessels +=1
    start_node = len(nodes)
    with open(filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader: 
            if line_count > 0:
                nodes.append((float(row[1])-FOV/2., float(row[2])-FOV/2., 0.000)) # x,y,z
            line_count +=1
        NSegments +=(int(row[0])-1)
    end_node = len(nodes)-1
    Vessels.append((start_node, end_node))

# Write connectivity data
for i in range(NVessels):
    for segment in range(Vessels[i][0], Vessels[i][1]-1):
        connectivity.append([segment, segment+1]) # Each node is linked to the next (no branching)
    connectivity.append([Vessels[i][1]-1, Vessels[i][1]]) # Link the terminal node to the vessel
    

# Create a root tree to link the vessels to a single inlet
root_nodes = [(-2.0*FOV, 0.000, 0.000), (-1.5*FOV, 0.000, 0.000)] # Root segment
root_connectivity = [[NSegments, -1]]                             # Connectivity for root tree
NSegments+=1

# Add a number of nodes along all sides of the square to connect the segmented vessels to the root

        



print('Vessels info:', Vessels)
print('Number of vessels:', NVessels)
print('Number of segments:', NSegments)
print('Number of nodes:', len(nodes))
print('Connectivity:', connectivity)
    
