'''
Tools for analysis of a .cco file
'''
import matplotlib.pyplot as plt
import numpy as np

def ReadCCO(ccoFileName):
    with open(ccoFileName, 'r') as f:
        # Unused *Tree information
        row = f.readline()
        row = f.readline()
        row = f.readline()

        # Vessels information
        row = f.readline()
        # print('Reading', row.strip(), '...')
        nVessels = int(f.readline())
        print("Reading", ccoFileName, "with", nVessels, "vessels in the tree (including root).")
        # print(nVessels, "vessels in the tree.")

        vessels = {}
        for i in range(nVessels):
            row = (f.readline()).split() # Split all columns in a list
            Id, xProx, xDist, r, q, stage = int(row[0]), [float(xi) for xi in row[1:4]], [float(xi) for xi in row[4:7]], float(row[12]), float(row[10]), int(row[-1])
            l = sum([(a-b)**2 for a,b in zip(xProx, xDist)])**.5

            vessels[Id] = [r, l, q, stage]

        row = f.readline()
        row = f.readline()
        # print('Reading', row.strip())
        connectivity = {}
        for i in range(nVessels):
            row = (f.readline().split())
            if (len(row)==2):
                Id, parentId = int(row[0]), int(row[1])
                connectivity[Id] = [parentId, []]
            else:
                Id, parentId, children = int(row[0]), int(row[1]), [int(i) for i in row[2:]]
                connectivity[Id] = [parentId, children]

        return vessels, connectivity

def VesselDensity(vessels, domainAreaInCm=0.09, upToStage=100):
    '''
    domainAreaInCm is the surface of the perfusion domain, i.e. the square 
    including the FAZ.
    upToStage is used to compute vessel density of the tree
    starting up to certain stage (by default includes all stages).
    Ex: VesselDensity(vessels, domainAreaInCm, 2) will return the 
    vessel density of the tree, excluding vessels that are grown 
    at stage 3 or later.
    '''

    VAD = 0.0
    vesselCount   = 0
    for Id in vessels:
        r, l, stage = vessels[Id][0], vessels[Id][1], vessels[Id][-1]
        if stage < upToStage and stage>=-1:
            VAD += 2.00*r*l # Surface of the vessel on a 2d projection
            vesselCount +=1
    VAD/=domainAreaInCm
    return VAD, vesselCount

def VesselPerimeterIndex(vessels, domainAreaInCm=0.09, upToStage=100):
    '''
    domainAreaInCm is the surface of the perfusion domain, i.e. the square 
    including the FAZ.
    upToStage is used to compute vessel density of the tree
    starting up to certain stage (by default includes all stages).
    Ex: VesselPerimeterIndex(vessels, domainAreaInCm, 2) will return the 
    vessel density of the tree, excluding vessels that are grown 
    at stage 3 or later.
    '''
    
    VPI = 0.0
    vesselCount   = 0
    for Id in vessels:
        r, l, stage = vessels[Id][0], vessels[Id][1], vessels[Id][-1]
        if stage < upToStage and stage>=-1:
            VPI += 2.00*l # Perimeter of the vessel segment on a 2D projection
            vesselCount +=1
    VPI/=domainAreaInCm
    VPI/=1e4                    # Convert from cm-1 to micron-1
    return VPI, vesselCount

def VesselDiameterIndex(vessels, domainAreaInCm=0.09, upToStage=100):
    '''
    domainAreaInCm is the surface of the perfusion domain, i.e. the square 
    including the FAZ.
    upToStage is used to compute vessel density of the tree
    starting up to certain stage (by default includes all stages).
    Ex: VesselDiameterIndex(vessels, domainAreaInCm, 2) will return the 
    vessel density of the tree, excluding vessels that are grown 
    at stage 3 or later.
    '''

    totalLength, totalArea = 0.0, 0.0
    vesselCount   = 0
    for Id in vessels:
        r, l, stage = vessels[Id][0], vessels[Id][1], vessels[Id][-1]
        if stage < upToStage and stage>=-1:
            totalLength += l
            totalArea   += l*2*r
            vesselCount +=1
    VDI = totalArea/totalLength*1e4 # Convert from cm to micron
    return VDI, vesselCount

def VesselComplexityIndex(vessels, domainAreaInCm=0.09, upToStage=100):
    '''
    domainAreaInCm is the surface of the perfusion domain, i.e. the square 
    including the FAZ.
    upToStage is used to compute vessel density of the tree
    starting up to certain stage (by default includes all stages).
    Ex: VesselComplexityIndex(vessels, domainAreaInCm, 2) will return the 
    vessel density of the tree, excluding vessels that are grown 
    at stage 3 or later.
    '''

    totalPerimeter, totalArea = 0.0, 0.0
    vesselCount   = 0
    for Id in vessels:
        r, l, stage = vessels[Id][0], vessels[Id][1], vessels[Id][-1]
        if stage < upToStage and stage>=-1:
            totalPerimeter += l*2
            totalArea      += l*2*r
            vesselCount    +=1
    VCI = (totalPerimeter**2)/(4*np.pi*totalArea)
    return VCI, vesselCount

def StatisticsMultipleTrees(ccoFiles, domainAreaInCm=0.09, upToStage=100):
    '''
    Compute vessel density mean and std 
    from all the trees in ccoFiles.
    '''

    VAD, VPI, VDI, VCI = [], [], [], []
    for ccoFile in ccoFiles:
        vessels, _ = ReadCCO(ccoFile)
        vesselDensity, vesselCount = VesselDensity(vessels, domainAreaInCm, upToStage)
        VAD.append(vesselDensity)            
        VPI.append(VesselPerimeterIndex(vessels, domainAreaInCm, upToStage)[0])
        VDI.append(VesselDiameterIndex(vessels, domainAreaInCm, upToStage)[0])
        VCI.append(VesselComplexityIndex(vessels, domainAreaInCm, upToStage)[0])
        print(vesselCount, "vessels included in the count.")

    return np.array(VAD), np.array(VPI), np.array(VDI), np.array(VCI)
        
        
