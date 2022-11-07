# %%
import pandas as pd
import numpy as np
from math import acos
from scipy.linalg import solve


# %%
def VerticalProjection(P):
    r = 1.14e04 # In microns
    Q = P
    Q[1:,-1]  = -(r**2 - P[1:,0]**2 - P[1:,1]**2)**0.5
    Q[0,:-1]  -= r  
    return Q

def StereographicProjection(P):
    r = 1.14e04
    Q = np.zeros_like(P)
    for i in range(P.shape[0]):
        if P[i,2] < 0:
            Q[i,:]  = [P[0,0], P[0,1], P[0,2]-r]
        else:
            h = P[i,0]**2 + P[i,1]**2
            Q[i,0] = r* 2*P[i,0]*r/(r**2+h)
            Q[i,1] = r* 2*P[i,1]*r/(r**2+h)
            Q[i,2] = r* (h-r**2)/(h+r**2)
    return Q
    
# %%
def Viscosity(d, hd=0.45, Type='Pries'):
    '''
    Diameter (d) should be in microns
    '''
    if Type=='Constant':
        return 3.6 # cP
    if Type=='Haynes':
        # Both muInf and delta are taken from Takahashi's model
        muInf = 1.09*np.exp(0.024*hd)
        delta = 4.29
        return muInf/( (1+delta/(d/2.0))**2 )
    else:
        mu045 = 220*np.exp(-1.3*d) + 3.2 - 2.44*np.exp(-0.06*d*0.645)
        C = (0.8 + np.exp(-0.075*d))*(-1+(1+10**-11*(d)**12)**-1)+(1+10**-11*(d)**12)**-1
        return ( 1 + (mu045-1)*((1-hd)**C-1)/((1-0.45)**C-1) )



# %%
def SolveFlow(Connectivity, Radii, Length, D, df, viscosityModel='Pries', Qperf=15.0):
    
    muLmin_to_m3sec = 1e-9/60
    muLmin_to_mum3sec = 1e9/60
    mmHg_to_Pa = 133.3224 # To Pa=kg.m^-1.s^-2
    cP_to_Pas  = 1e-3     # To Pa.s=kg.m^-1.s^-1
    micron_to_meter = 1e-6
    mum_to_mm = 1e-3
    Pa_to_dynecm2 = 1e1

    narcs, nnodes = Connectivity.shape
    R = np.zeros((narcs,narcs))
    for i in range(narcs):
        r,l = Radii[i], Length[i]  # In meter
        mu  = Viscosity(2*r, Type=viscosityModel) * cP_to_Pas 
        # print(f'{r=} and {mu=}')
        R[i,i] = 8.0 * mu * l*micron_to_meter/(np.pi * ((r*micron_to_meter)**4))
    
    # Add the boundary conditions
    qperf = Qperf * muLmin_to_m3sec # Qperf is in muL/min, changed to m^3/s
    qcap  = qperf / df.value_counts('Boundary')[-1]

    # RHS for flow rates: qbar = qperf if inlet node, qcap otherwise
    qbar = np.ones((nnodes,)) # Flow rate at nodes
    qbar = np.matmul(D,qbar)
    qbar[qbar>0] = qperf 
    qbar[qbar<0] = qcap

    pbar = np.ones((nnodes,))
    pbar = D.dot(pbar)
    pbar[pbar>0] = 50 * mmHg_to_Pa
    pbar[pbar<0] = 25 * mmHg_to_Pa
    
    D = np.abs(D)
    # D[-1,-1] = 0 # Setting this entry to 0 allows to specify flow at inlet 
    I = np.identity(D.shape[0])
    
    M = np.block([[
    R, -Connectivity],
    [(I-D).dot(Connectivity.T), D]])
    
    RHS = np.concatenate([np.zeros(narcs,), (I-D).dot(qbar)+D.dot(pbar)])
    x = solve(M,RHS)
    f,p = np.abs(x[:narcs]), x[narcs:] / mmHg_to_Pa
    dp  = np.abs(Connectivity.dot(p)) / mmHg_to_Pa   
    
    for i in range(f.shape[0]):
        df.at[i,'CompFlow'] = df.at[i, 'Flow']
        df.at[i,'Flow'] = f[i] / muLmin_to_m3sec
        df.at[i,'PressureDrop'] = dp[i] 
        df.at[i,'Pressure'] = p[df.at[i,'DistalNode']]  + dp[i] 
    
    df['Viscosity'] = Viscosity(df.Radius.to_numpy()*2, Type=viscosityModel) # In cP
    df['FlowVelocity'] = df.Flow.to_numpy() * muLmin_to_mum3sec/(np.pi*(df.Radius.to_numpy()**2)) * mum_to_mm # Converts flow in muL/min to mum^3/min then mm/s
    df['WSR'] = 4*df.Flow.to_numpy() * muLmin_to_mum3sec /(np.pi*df.Radius.to_numpy()**3)  
    df['WSS'] = df.Viscosity.to_numpy()*df.WSR.to_numpy() * cP_to_Pas * Pa_to_dynecm2
    df['Diameter'] = df.Radius.to_numpy()*2.0
    df['Resistance'] = 8.0*df.Viscosity.to_numpy()*df.Length.to_numpy()/(np.pi*(df.Radius.to_numpy()**4)) * cP_to_Pas / (micron_to_meter**3) # Results in Pa.s.m^-3
    # df['Resistance'] = df['Pressure'].max()/df['Flow'] 
    
    cols = ['Flow', 'Pressure', 'FlowVelocity', 'PressureDrop',
            'Boundary', 'Viscosity', 'WSR', 'WSS', 'Radius', 'Length', 
            'Bifurcation', 'Stage', 'Id',
            'ParentId', 'BranchesId',
            'xProx', 'xDist',
            'DistalNode', 'Diameter', 'Resistance']
    dfNew = df[cols]
    dfNew.set_index('Id', inplace=True, drop=False)
    
    print(f'Flow loss (C*f).sum() = {Connectivity.T.dot(df.Flow.to_numpy()).sum()}, with q_in={dfNew.Flow.max()}')
    print(f"Flow loss q_in-q_out = {df.loc[df.Boundary==1,'Flow'].sum() - df.loc[df.Boundary==-1, 'Flow'].sum()}")
    
    
    return dfNew

# %%
def ReadCCO(CCOFile, project=True, treeType='Object'):
    
    SegmentsData = []
    DistalNodes  = dict()
    
    nodeName = -1

    with open(CCOFile, 'r') as f:
            
        print(f.readline().strip())
        token = f.readline().split()
        inletNode = ([float(xi) for xi in token[:3]])
        print(f"Inlet at {inletNode}")

        f.readline()
        print(f.readline().strip())
        nVessels = int(f.readline())
        print(f'The tree has {nVessels} vessels.')
        
        for i in range(nVessels):
            
            row = (f.readline()).split()
            
            nodeName+=1
            DistalNodes[int(row[0])] = nodeName
                
            if treeType=='Object':
                # Uncomment to keep vessels added in specific stages
                if True:
                # if int(row[-1]) > -2: 

                    # Id, xProx, xDist, radius, length, flow computed by CCO, distalNode, stage
                    SegmentsData.append([int(row[0]),
                    np.array([float(x) for x in row[1:4]])*1e4, 
                    np.array([float(x) for x in row[4:7]])*1e4,
                    float(row[12])*1e4,
                    np.linalg.norm(np.array([float(x) for x in row[1:4]])*1e4
                                -np.array([float(x) for x in row[4:7]])*1e4),
                    float(row[13]),
                    nodeName,
                    int(row[-1])])
                    

            else:
                # Id, xProx, xDist, radius, length, flow computed by CCO, distalNode, vesselId (for multiple segments vessels?)
                SegmentsData.append([int(row[0]),
                np.array([float(x) for x in row[1:4]])*1e4, 
                np.array([float(x) for x in row[4:7]])*1e4,
                float(row[8])*1e4,
                float(row[10])*1e4,
                float(row[11]),
                nodeName,
                int(row[-2])])
                
        if treeType=='Object':            
            df = pd.DataFrame(SegmentsData, columns=['Id', 'xProx', 'xDist', 'Radius', 'Length', 'Flow','DistalNode','Stage'])
        else:
            df = pd.DataFrame(SegmentsData, columns=['Id', 'xProx', 'xDist', 'Radius', 'Length', 'Flow','DistalNode','Stage'])

        df['Inlet'] = False
        df['Outlet'] = False
        df = df.set_index('Id', drop=False)
        df['ParentId'] = -1
        df['BranchesId'] = [[] for i in range(df.shape[0])]
    
        f.readline()
        print(f.readline().strip())
        
        NodesConnections = []
        SegNewName = -1
        for i in range(nVessels):
            row = (f.readline()).split()
            SegmentId, ParentId, BranchesIds = int(row[0]), int(row[1]), [int(x) for x in row[2:]]

            if SegmentId in DistalNodes:    
                ProximalNode = DistalNodes[SegmentId]
                branchesId = []
                for connection in BranchesIds:
                    DistalNode = DistalNodes[connection]
                    SegNewName +=1
                    NodesConnections.append((ProximalNode, DistalNode, SegNewName, connection))
                    branchesId.append(connection)
                df.at[SegmentId, 'BranchesId'] = branchesId    
                
                if not BranchesIds:
                    df.at[SegmentId, 'Outlet'] = True
                
                if not ParentId in DistalNodes: # Inlet node, need to add the proximal node to the tree
                    df.at[SegmentId, 'Inlet'] = True
                    nodeName+=1
                    SegNewName +=1
                    NodesConnections.append((nodeName, DistalNodes[SegmentId], SegNewName, SegmentId))
                else:
                    df.at[SegmentId, 'ParentId'] = ParentId

    
    ## Gives each inlet and its downstream branches a number
    def AssignBranchNumberToDownstreamVessels(InletIdx, branchNumber):
        
        # branchNumber = df.at[InletIdx, 'BranchNumber']
        branches = df.at[InletIdx, 'BranchesId']
        
        for branchId in branches:
            branchIdx = df[df.Id==branchId].index[0]
            df.at[branchIdx, 'BranchNumber'] = branchNumber
            df.at[branchIdx, 'Bifurcation']  = df.at[InletIdx,'Bifurcation']+1 
            AssignBranchNumberToDownstreamVessels(branchIdx, branchNumber)       
        
    df['BranchNumber'] = -1
    df['Bifurcation']  = 0
    for i,row in enumerate(df[df.Inlet].iterrows()):
        df.at[row[0], 'BranchNumber'] = i        
        AssignBranchNumberToDownstreamVessels(row[0],i)
                    
    # Project from the plane to the sphere
    if project:
        radiusEyeball = 23e3
        # xProx = ProjectOnSphere(np.vstack(df['xProx'].to_numpy().ravel()), r=radiusEyeball)
        # xDist = ProjectOnSphere(np.vstack(df['xDist'].to_numpy().ravel()), r=radiusEyeball)
        xProx = StereographicProjection(np.vstack(df.xProx.to_numpy().ravel()))
        xDist = StereographicProjection(np.vstack(df.xDist.to_numpy().ravel()))
        df['xProx'] = [xProx[i,:] for i in range(xProx.shape[0])] 
        df['xDist'] = [xDist[i,:] for i in range(xDist.shape[0])]
        df['Length'] = np.linalg.norm(xProx-xDist, axis=1)

    ## Create the matrices for the solver
    ConnectivityMatrix = np.zeros((nodeName+1, len(SegmentsData)))
    Radii  = np.zeros((len(SegmentsData),))
    Length = np.zeros((len(SegmentsData),))
    df['SegName'] = df.index
    df['Id'] = df.index
    df['Boundary'] = 0
    D = np.zeros((nodeName+1,nodeName+1)) # Decision matrix

    for proxNode, distNode, SegmentName, SegmentId in NodesConnections:
        
        ConnectivityMatrix[proxNode, SegmentName] = 1
        ConnectivityMatrix[distNode, SegmentName] = -1
        Radii[SegmentName] = df.at[SegmentId, 'Radius']
        xProx, xDist = df.at[SegmentId, 'xProx'], df.at[SegmentId, 'xDist']
        Length[SegmentName] = np.linalg.norm(xProx-xDist) 
        
        df.at[SegmentId, 'SegName'] = SegmentName
        if df.at[SegmentId, 'Inlet']:
            df.at[SegmentId,'Boundary'] = 1
            D[proxNode, proxNode] = 1
        elif df.at[SegmentId, 'Outlet']:
            df.at[SegmentId,'Boundary'] = -1
            D[distNode, distNode] = -1
    
    df = df.set_index('SegName')
            
    return ConnectivityMatrix.T, Radii, Length, D, df

