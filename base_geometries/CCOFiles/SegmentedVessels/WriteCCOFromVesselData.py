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
    if i>0:
        connectivity.append([Vessels[i][0], Vessels[i][0]+1])
    for segment in range(Vessels[i][0]+1, Vessels[i][1]):
        connectivity.append([segment, segment-1, segment+1]) # Each node is linked to the next (no branching)
    connectivity.append([Vessels[i][1], Vessels[i][1]-1])    # Link the terminal node to the vessel

    
# Create a root tree to link the vessels to a single inlet
root_nodes = [(-2.0*FOV, 0.000, 0.000), (-0.75*FOV, 0.000, 0.000)] # Root segment
root_connectivity = [[NSegments, -1]]                              # Connectivity for root tree
NSegments+=1

## Define the points that will be linked to the segmented vessels
n = 6                           # Number of nodes along each sides of the square
h = 1.3*FOV/n
x,y = root_nodes[1][:-1]        # Point of the bifurcation

### The branch going above the square
nNodesAboveBranch = 0
nn = int(n/2)
for i in range(nn):
    y+=h
    root_nodes.append((x,y,0.000))
    nNodesAboveBranch+=1
for i in range(n):
    x+=h
    root_nodes.append((x,y,0.000))
    nNodesAboveBranch+=1
for i in range(nn):
    y-=h
    root_nodes.append((x,y, 0.000))
    nNodesAboveBranch+=1

root_connectivity[0].append(NSegments)
for node in root_nodes[2:-1]:
    root_connectivity.append([NSegments, NSegments-1, NSegments+1])
    NSegments+=1
root_connectivity.append([NSegments, NSegments-1]) # Terminal vessel
NSegments+=1
    
### The branch going under the square
x,y = root_nodes[1][:-1]
nn = n-nn
for i in range(nn):
    y-=h
    root_nodes.append((x,y,0.000))
for i in range(n):
    x+=h
    root_nodes.append((x,y,0.000))
for i in range(nn-1):
    y+=h
    root_nodes.append((x,y, 0.000))
    
root_connectivity[0].append(NSegments)
root_connectivity.append([NSegments, root_connectivity[0][0], NSegments+1])
NSegments+=1
for node in root_nodes[2+nNodesAboveBranch:-1]:
    root_connectivity.append([NSegments, NSegments-1, NSegments+1])
    NSegments+=1
root_connectivity.append([NSegments, NSegments-1]) # Terminal vessel
NSegments+=1


def closestNode(IdClosestNodeInRootTree):
    # Needs to return the closest neigbhours to which connect the end points of both tree (root and segmented)
    # while not assigning a vessel from the root tree to more than one segmented vessel

def savetree(root_connectivity, root_nodes, connectivity, nodes):
    
    with open('RootTreeFromImage.cco', 'w') as f:
        # global tree information
        xPerf = root_nodes[0]
        qProx = 10.00
        psiFactor = 6.25e-07
        dp        = 30.4e+07
        nTerms    = 2 + NVessels
        refPressure = 9.80
        pointCounter = 0
        rootRadius   = Radius
        variationTolerance = 1e-05

        f.write('*Tree\n')
        f.write(f'{xPerf[0]} {xPerf[1]} {xPerf[2]} {qProx} {psiFactor} {dp} {nTerms} {refPressure} {pointCounter} {rootRadius} {variationTolerance}\n')
        
        vtkSegmentId   = 0
        xProx          = []
        xDist          = []
        unused         = 0.0
        qReservedFrac  = 0.0
        branchingMode  = 2      # BRANCHING_MODE {NO_BRANCHING, RIGID_PARENT, DEFORMABLE_PARENT, DISTAL_BRANCHING, ONLY_AT_PARENT_HOTSPOTS}
        vesselFunction = 1      # VESSEL_FUNCTION { DISTRIBUTION, PERFORATOR, TRANSPORT }
        stage          = 1      # 0 for the 'root', 1 for the segmented vessels

        # Writing vessels data 
        f.write('\n*Vessels\n')
        f.write(str(len(root_connectivity) + len(connectivity)) + '\n')
        xProx = nodes[0]
        for node in nodes[1:]:
            xDist = node
            f.write(f'{vtkSegmentId:d} {xProx[0]} {xProx[1]} {xProx[2]} {xDist[0]} {xDist[1]} {xDist[2]} {unused} {unused} {unused} {qReservedFrac} {branchingMode:d} {Radius} {unused} {unused} {unused} {unused} {vesselFunction} {unused} {unused} {stage:d}\n')
            xProx = xDist
            vtkSegmentId+=1

        xProx = root_nodes[0]
        stage = 0
        for node in root_nodes[1:]:
            xDist = node
            f.write(f'{vtkSegmentId:d} {xProx[0]} {xProx[1]} {xProx[2]} {xDist[0]} {xDist[1]} {xDist[2]} {unused} {unused} {unused} {qReservedFrac} {branchingMode:d} {Radius} {unused} {unused} {unused} {unused} {vesselFunction} {unused} {unused} {stage:d}\n')
            xProx = xDist
            vtkSegmentId+=1

        # Identify the points to connect the root to the segmented vessels
        IdStartSegementedVessel = []
        IdClosestNodeInRootTree = []
        for vessel in Vessels:
            IdStartSegmentedVessel.append(vessel[0]) # Id of first segment of each segmented vessel
            IdClosestNodeInRootTree.append(closestNode(IdClosestNodeInRootTree)) # Which node to link it too within the root tree
        
        # Add the connecting vessels and update the connectivity table
        print('Length of IdStartSegmentedVessel, IdClosestNodeInRootTree, should be equal:', len(IdStartSegmentedVessel), len(IdClosestNodeInRootTree))
        for IdProx, IdDist in zip(IdClosestNodeInRootTree, IdStartSegmentedVessel):
            xProx = root_nodes[IdProx]
            xDist = nodes[IdDist]

            root_connectivity.append([vtkSegmentId, IdProx-1, IdDist-1])
            connectivity[IdDist-1].insert(1, IdProx-1)

            f.write(f'{vtkSegmentId:d} {xProx[0]} {xProx[1]} {xProx[2]} {xDist[0]} {xDist[1]} {xDist[2]} {unused} {unused} {unused} {qReservedFrac} {branchingMode:d} {Radius} {unused} {unused} {unused} {unused} {vesselFunction} {unused} {unused} {stage:d}\n')
            vtkSegmentId+=1

        # Write the connectivity table
        f.write('\n*Connectivity\n')
        for vessel in connectivity:
            string = ' '.join(str(Id) for Id in vessel)
            f.write(string + '\n')
        for vessel in root_connectivity:
            string = ' '.join(str(Id) for Id in vessel)
            f.write(string + '\n')


savetree(root_connectivity, root_nodes, connectivity, nodes)

