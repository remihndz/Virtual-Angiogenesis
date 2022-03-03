'''
Tools for the computation of 
Horton-Strahler stream order.
'''

import matplotlib.pyplot as plt
import numpy as np
import copy

def ReadTree(ccoFile):
    with open(ccoFile, 'r') as f:
        # Unused '*Tree' information
        row = f.readline()
        row = f.readline()
        row = f.readline()

        # Vessels data
        print('Reading', f.readline().strip(), '...')
        nVessels = int(f.readline())
        print("Number of vessels in the tree: ", nVessels)

        treeData = []
        for i in range(nVessels):
            row = (f.readline()).split() # Split all columns in a list
            Id, xProx, xDist, r, q, stage = int(row[0]), row[1:4], row[4:7], float(row[12]), float(row[10]), int(row[-1])
            l = sum([(float(a)-float(b))**2 for a,b in zip(xProx, xDist)])**.5
            treeData.append([Id, r, l, q, stage])
        
        row = f.readline()
        print('Reading', f.readline().strip())
        treeConnectivity = []
        for i in range(nVessels):
            row = (f.readline()).split()
            if (len(row)==2):
                Id, parentId = int(row[0]), int(row[1])
                treeConnectivity.append([Id, parentId, []])
            else:
                Id, parentId, children = int(row[0]), int(row[1]), [int(i) for i in row[2:]]
                treeConnectivity.append([Id, parentId, children])         

    return treeData, treeConnectivity

def StrahlerOrder(treeData, treeConnectivity):

    treeConnectivity.sort()
    currentOrder = 0
    vessel = treeConnectivity[0]
    vesselDict = {vessel[0]: [vessel[1], vessel[2], currentOrder]}
    for vessel in treeConnectivity[1:]:
        if len(vessel)>2:
            IdVessel, IdParent, ChildrenList = vessel
        else:
            IdVessel, IdParent = vessel
            ChildrenList = []
        vesselDict[IdVessel] = [IdParent, ChildrenList, currentOrder]
    prunedTree = copy.deepcopy(vesselDict)
    
    def CurrentLeaves(prunedTree):
        leavesList = []
        for Id, vessel in prunedTree.items():
            if vessel[-1] == 0 and len(vessel[-2])==0:
                leavesList.append(Id)
        return leavesList

    def UpdateAndPrune(prunedTree, leaves):        
        # Assign an order to the current leaves
        copyTree = copy.deepcopy(prunedTree)
        increment= False

        for leaf in leaves:
            Children = treeConnectivity[leaf][-1]
            if len(Children)==2 and copyTree[Children[0]][-1]==currentOrder and copyTree[Children[1]][-1]==currentOrder:
                copyTree[leaf][-1] = currentOrder+1
                increment = True
                
            else:
                childrenOrders = [0]
                for child in Children:
                    childrenOrders.append(copyTree[child][-1])
                ''' 
                Two different rules, the max rule is the one from Wikipedia and gives fewer orders.
                Both give the right results for the test tree.
                '''
                copyTree[leaf][-1] = max(1,max(childrenOrders)) # currentOrder

        # Prune
        for leaf in leaves:
            parentId = copyTree[leaf][0]
            if parentId>=0:
                copyTree[parentId][1].remove(leaf)
        return copyTree, increment

    currentOrder = 1
    computeForTheRoot = iter([True, False]) # Just used to run one more time once the root is reached
    while prunedTree[0][-1]==0 or next(computeForTheRoot):
        leafList = CurrentLeaves(prunedTree)
        prunedTree, increment = UpdateAndPrune(prunedTree, leafList)
        if increment: currentOrder+=1
        
        if 0 in leafList:
            break

    # Add stream order to vessel data. Get rid of the root vessels that are outside the domain
    finalTree = []
    for vessel in treeData:
        stage = vessel[-1]
        if stage >= -1:
            Id = int(vessel[0])
            order = prunedTree[Id][-1]
            treeData[Id].append(order)
            finalTree.append(treeData[Id])

    # Shift the indexes to compensate for the vessels that were deleted
    nVesselsDeleted = len(treeData)-len(finalTree)
    for i in range(len(finalTree)):
        finalTree[i][0] -= nVesselsDeleted
        
    return finalTree

def TreeStatistics(orderedTree, treeData):
    return #stats

def PlotTreeStatistics(orderedTree, outputImageFile=None):

    maxOrder = np.array(orderedTree)[:,-1].astype(int).max()
    data = np.array(orderedTree)

    dataRadius = np.zeros((maxOrder, len(orderedTree)))
    dataLength = np.zeros((maxOrder, len(orderedTree)))
    dataVolume   = np.zeros((maxOrder, len(orderedTree)))
    dataAspectRatio   = np.zeros((maxOrder, len(orderedTree)))
    
    for Id, radius, length, flow, stage, order in orderedTree:
        r,l = radius*1e4, length*1e4 # Convert to microns
        i,j = int(order)-1, int(Id)
        dataRadius[i, j] = r
        dataLength[i, j] = l
        dataVolume[i, j] = r*r*l*np.pi
        dataAspectRatio[i, j] = l/r
        
    meanRadius = np.mean(dataRadius, axis = 1, where=dataRadius>0)
    stdRadius  = np.std(dataRadius, axis = 1, where=dataRadius>0)
    meanLength = np.mean(dataLength, axis = 1, where=dataRadius>0)
    stdLength  = np.std(dataLength, axis = 1, where=dataRadius>0)
    meanVolume = np.mean(dataVolume, axis = 1, where=dataRadius>0)
    stdVolume  = np.std(dataVolume, axis = 1, where=dataRadius>0)
    meanAspectRatio = np.mean(dataAspectRatio, axis = 1, where=dataRadius>0)
    stdAspectRatio  = np.std(dataAspectRatio, axis = 1, where=dataRadius>0)
    
    x = np.flip(np.arange(maxOrder) + 1)
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col')
    fig.set_figheight(7.2)
    fig.set_figwidth(8.8)

    scale = 'linear' # 'log'
    ax1.errorbar(x, meanRadius, stdRadius, marker='s', fillstyle='full')
    ax1.set_yscale(scale)
    ax1.set_title('Radius vs. Strahler order')
    ax2.errorbar(x, meanLength, stdLength, marker='s')
    ax2.set_yscale(scale)
    ax2.set_title('Length vs. Strahler order')
    ax3.errorbar(x, meanVolume, stdVolume, marker='s')
    ax3.set_yscale(scale)
    ax3.set_title('Volume vs. Strahler order')
    ax4.errorbar(x, meanAspectRatio, stdAspectRatio, marker='s')
    ax4.set_yscale(scale)
    ax4.set_title('Aspect ratio vs. Strahler order')

    fig2, axes = plt.subplots(2,2, sharex='col')
    fig2.set_figheight(7.2)
    fig2.set_figwidth(8.8)
    scale = 'linear'
    axes[0,0].boxplot(dataRadius.transpose())
    axes[0,0].set_yscale(scale)
    axes[0,0].set_title('Radius vs. Strahler order')
    axes[0,1].boxplot(dataLength.transpose())
    axes[0,1].set_yscale(scale)
    axes[0,1].set_title('Length vs. Strahler order')
    axes[1,0].boxplot(dataVolume.transpose())
    axes[1,0].set_yscale(scale)
    axes[1,0].set_title('Volume vs. Strahler order')
    axes[1,1].boxplot(dataAspectRatio.transpose())
    axes[1,1].set_yscale(scale)
    axes[1,1].set_title('Aspect ratio vs. Strahler order')


    if outputImageFile:
        np.savetxt(outputImageFile[:-4] + 'vessels.dat', data)
        print('\nWriting outputs (image+data) in', outputImageFile[:-4] + '.png/.dat')        
        data = np.concatenate((x.reshape((-1,1)), meanRadius.reshape((-1,1)), stdRadius.reshape((-1,1)), meanLength.reshape((-1,1)), stdLength.reshape((-1,1)), meanVolume.reshape((-1,1)), stdVolume.reshape((-1,1)), meanAspectRatio.reshape((-1,1)), stdAspectRatio.reshape((-1,1))), axis=1)
        #np.savetxt(outputImageFile[:-4] + '.dat', data)
        #plt.savefig(outputImageFile)
    plt.show()
    
    return

