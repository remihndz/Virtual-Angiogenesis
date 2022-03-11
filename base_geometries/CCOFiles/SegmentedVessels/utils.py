from operator import itemgetter

def CreateRootTree(FOV, Radius, n):

    root_nodes = [(-2.0*FOV, 0.000, 0.000), (-0.75*FOV, 0.000, 0.000)] # Root segment
    RootConn = []                              # Connectivity for root tree
    RootData = []

    ## Define the points that will be linked to the segmented vessels
    n = 8                           # Number of nodes along each sides of the square
    h = 1.3*FOV/n
    x,y = root_nodes[1][:-1]        # Point of the bifurcation
    
    ### The branch going above the square
    nNodesAboveBranch = 0
    nn = int(n/2)
    for i in range(nn):
        y+=h
        root_nodes.append((x,y,0.000))
        nNodesAboveBranch+=1
    for i in range(n+1):
        x+=h
        root_nodes.append((x,y,0.000))
        nNodesAboveBranch+=1
    for i in range(nn):
        y-=h
        root_nodes.append((x,y, 0.000))
        nNodesAboveBranch+=1
        
    vtkSegmentId = 0
    for xProx, xDist in zip(root_nodes[:-2], root_nodes[1:-1]):

        RootData.append([vtkSegmentId, xProx[0], xProx[1], xProx[2], xDist[0], xDist[1], xDist[2]])
        RootConn.append([vtkSegmentId, vtkSegmentId-1, vtkSegmentId+1])
        vtkSegmentId+=1

    xProx, xDist = root_nodes[-2], root_nodes[-1]
    RootData.append([vtkSegmentId, xProx[0], xProx[1], xProx[2], xDist[0], xDist[1], xDist[2]])
    RootConn.append([vtkSegmentId, vtkSegmentId-1])
    vtkSegmentId+=1
    RootConn[0]  = [0, -1, 1, vtkSegmentId]
    RootConn[1]  = [1, 0, 2]

    ### The branch going under the square
    x,y = root_nodes[1][:-1]
    root_nodes = [(x,y, 0.000)]

    nn = n-nn
    for i in range(nn):
        y-=h
        root_nodes.append((x,y,0.000))
    for i in range(n+1):
        x+=h
        root_nodes.append((x,y,0.000))
    for i in range(nn-1):
        y+=h
        root_nodes.append((x,y, 0.000))

    xProx, xDist = root_nodes[0], root_nodes[1]
    RootData.append([vtkSegmentId, xProx[0], xProx[1], xProx[2], xDist[0], xDist[1], xDist[2]])        
    RootConn.append([vtkSegmentId, RootConn[0][0], vtkSegmentId+1])
    vtkSegmentId+=1

    for xProx, xDist in zip(root_nodes[1:-1], root_nodes[2:]):
        RootData.append([vtkSegmentId, xProx[0], xProx[1], xProx[2], xDist[0], xDist[1], xDist[2]])        
        RootConn.append([vtkSegmentId, vtkSegmentId-1, vtkSegmentId+1])
        vtkSegmentId+=1

    RootConn[-1] = [vtkSegmentId-1, vtkSegmentId-2]
    return RootData, RootConn, vtkSegmentId


def LinkToRootTree(VesselsData, VesselsConn, RootData, RootConn):
    # Needs to return the closest neigbhours to which connect the end points of both tree (root and segmented)
    # while not assigning a vessel from the root tree to more than one segmented vessel
    def dist(xDist, xProx):
        return sum([(x-y)**2 for x,y in zip(xDist, xProx)])**.5

    links = []
    UsedSegments = [0]
    
    for i,vessel_segment in enumerate(VesselsData):
        SegmentId = vessel_segment[0]
        if VesselsConn[i][1] == -2: # Identifies a vessel that needs to be connected

            xDist = vessel_segment[4:7]
            distToxDist = []
            
            for segment in RootData:
                if not segment[0] in UsedSegments:
                    xProx = segment[4:7]
                    distToxDist.append((dist(xDist, xProx), segment[0]))

            distToxDist = sorted(distToxDist, key=itemgetter(0))
            IdClosestNeighbour = distToxDist[0][1]
            links.append([i, IdClosestNeighbour])
            UsedSegments.append(IdClosestNeighbour)

    # Update the connectivity tables and add the connecting segments in the data tables
    vtkSegmentId = VesselsData[-1][0]+1
    for IdVessel, IdRoot in links:
        ## Add the segments data to the root tree
        xProx = RootData[IdRoot][4:7]
        xDist = VesselsData[IdVessel][1:4]
        RootData.append([vtkSegmentId, xProx[0], xProx[1], xProx[2], xDist[0], xDist[1], xDist[2]])

        ## Update the connectivity tables
        VesselsConn[IdVessel][1] = vtkSegmentId
        RootConn[IdRoot].append(vtkSegmentId)
        RootConn.append([vtkSegmentId, IdRoot, VesselsConn[IdVessel][0]])
        vtkSegmentId +=1

    return    

def SaveTree(filename, VesselsData, VesselsConn, RootData, RootConn):
    Radius = 0.004
    # global tree information
    xPerf = RootData[0][1:4]
    qProx = 9.41e-6             # in cl/s from Takahashi
    refPressure = 22.9          # in mmHg from Takahashi
    psiFactor = 6.25e-07
    dp        = 3.04e+07
    nTerms    = 2
    for segment in VesselsData:
        if len(segment) == 2:
            nTerms+=1

    pointCounter = 0
    rootRadius   = Radius*2
    variationTolerance = 1e-05

    with open(filename, 'w') as f:
        f.write('*Tree\n')
        f.write(f'{xPerf[0]} {xPerf[1]} {xPerf[2]} {qProx} {psiFactor} {dp} {nTerms} {refPressure} {pointCounter} {rootRadius} {variationTolerance}\n\n')
        
        # Vessel information
        unused         = 0.0
        qReservedFrac  = 0.0
        branchingMode  = 0 
        vesselFunction = 2      # VESSEL_FUNCTION { DISTRIBUTION, PERFORATOR, TRANSPORT }
        stage          = -2     # -2 for the 'root', -1 for the segmented vessels
        
        f.write('*Vessels\n')
        f.write(str(len(VesselsData) + len(RootData)) + '\n')
        for segment in RootData:
            vtkSegmentId = segment[0]
            xProx = segment[1:4]
            xDist = segment[4:7]
            f.write(f'{vtkSegmentId:d} {xProx[0]} {xProx[1]} {xProx[2]} {xDist[0]} {xDist[1]} {xDist[2]} {unused} {unused} {unused} {qReservedFrac} {branchingMode:d} {Radius} {unused} {unused} {unused} {unused} {vesselFunction} {unused} {unused} {stage:d}\n')

        
        stage+=1
        branchingMode = 1       # BRANCHING_MODE {NO_BRANCHING, RIGID_PARENT, DEFORMABLE_PARENT, DISTAL_BRANCHING, ONLY_AT_PARENT_HOTSPOTS}
        vesselFunction = 1      # VESSEL_FUNCTION { DISTRIBUTION, PERFORATOR, TRANSPORT }
        for segment in VesselsData:
            vtkSegmentId = segment[0]
            xProx = segment[1:4]
            xDist = segment[4:7]
            
            f.write(f'{vtkSegmentId:d} {xProx[0]} {xProx[1]} {xProx[2]} {xDist[0]} {xDist[1]} {xDist[2]} {unused} {unused} {unused} {qReservedFrac} {branchingMode:d} {Radius} {unused} {unused} {unused} {unused} {vesselFunction} {unused} {unused} {stage:d}\n')

        f.write('\n*Connectivity\n')
        for segment in RootConn:
            s = ' '.join(str(Id) for Id in segment)
            f.write(s+'\n')
        for segment in VesselsConn:
            s = ' '.join(str(Id) for Id in segment)
            f.write(s+'\n')


