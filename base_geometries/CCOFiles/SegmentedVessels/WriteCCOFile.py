import sys
import glob
import csv
import os
from utils import CreateRootTree, LinkToRootTree, SaveTree

patientFolder = sys.argv[-1]
FOV    = 0.3                    # Size of the field of view. Shifts the FAZ to (x,y,z)=(0,0,0), in cm
Radius = 0.005                  # Fixed radius for all vessels, in cm

                                # TODO: add radius for each vessel
NVessels = 0                    # Number of vessels
NSegments = 0                   # Total number of segments in the segmented vessels
Vessels  = []                   # Container for vessels [(start_node, end_node)]

nodes = []                      # Container for nodes [(x0,y0,z0), (x1,y1,z1)...]
connectivity = []               # Container for connectivity 

## Save each segments data to a container. Each entry is a line in the cco file
VesselsData = []
VesselsConn = []

# Create a root tree
n = 6
RootData, RootConn, vtkSegmentId = CreateRootTree(FOV, Radius, n)

# Read data from files
for filename in glob.iglob(patientFolder + '/Vessel*'):
    NVessels +=1
    with open(filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader)
        row = next(csv_reader)
        xProx = (float(row[1])-FOV/2., float(row[2])-FOV/2., 0.000)
        isFirstSegment = True
        for row in csv_reader:
            xDist = (float(row[1])-FOV/2., float(row[2])-FOV/2., 0.000) # x,y,z

            ## Store vessel data
            VesselsData.append([vtkSegmentId, xProx[0], xProx[1], xProx[2], xDist[0], xDist[1], xDist[2]])

            ## Store vessel connectivity
            if isFirstSegment:
                VesselsConn.append([vtkSegmentId, -2]) # -2 identifies a segment that needs a parent in the root tree
                isFirstSegment=False
            else:
                VesselsConn[-1].append(vtkSegmentId)
                VesselsConn.append([vtkSegmentId, vtkSegmentId-1])
            vtkSegmentId+=1
            xProx = xDist

LinkToRootTree(VesselsData, VesselsConn, RootData, RootConn)
patientName = input("Name of the output cco file (without extension):")
SaveTree('../' + patientName+'.cco', VesselsData, VesselsConn, RootData, RootConn)
print('Tree save in: ' + '../' + patientName+'.cco')
