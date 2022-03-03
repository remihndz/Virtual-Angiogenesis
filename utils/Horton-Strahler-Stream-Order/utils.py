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
            treeData.append([Id, r*1e4, l*1e4, q, stage]) # Convert to mm
        
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
                copyTree[leaf][-1] = currentOrder # max(1,max(childrenOrders)) 

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

    # Add stream order to vessel data
    treeDataNoRoot = []
    count = 0
    for vessel in treeData:
        Id = int(vessel[0])
        order = prunedTree[Id][-1]
        treeData[Id].append(order)

        if treeData[Id][-2]>=-2:
            Id, radius, length, flow, stage, order = treeData[Id]
            treeDataNoRoot.append([count, radius, length, flow, order])
            count+=1
    return treeDataNoRoot

def TreeStatistics(orderedTree, treeData):
    return #stats

def PlotTreeStatistics(orderedTree, outputImageFile=None):

    maxOrder = np.array(orderedTree)[:,-1].astype(int).max()

    dataRadius = np.zeros((maxOrder, len(orderedTree)))
    dataLength = np.zeros((maxOrder, len(orderedTree)))
    dataVolume   = np.zeros((maxOrder, len(orderedTree)))
    dataAspectRatio   = np.zeros((maxOrder, len(orderedTree)))
    orderDistribution = np.zeros((maxOrder,))
    
    for Id, radius, length, flow, order in orderedTree:
        i,j = int(order)-1, int(Id)
        dataRadius[i, j] = radius
        dataLength[i, j] = length
        dataVolume[i, j] = radius*radius*length*np.pi
        dataAspectRatio[i, j] = length/radius
        orderDistribution[i] += 1
        
    meanRadius = np.mean(dataRadius, axis = 1, where=dataRadius>0)
    stdRadius  = np.std(dataRadius, axis = 1, where=dataRadius>0)
    meanLength = np.mean(dataLength, axis = 1, where=dataRadius>0)
    stdLength  = np.std(dataLength, axis = 1, where=dataRadius>0)
    meanVolume = np.mean(dataVolume, axis = 1, where=dataRadius>0)
    stdVolume  = np.std(dataVolume, axis = 1, where=dataRadius>0)
    meanAspectRatio = np.mean(dataAspectRatio, axis = 1, where=dataRadius>0)
    stdAspectRatio  = np.std(dataAspectRatio, axis = 1, where=dataRadius>0)

    x = np.arange(maxOrder) + 1 
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col')
    fig.set_figheight(7.2)
    fig.set_figwidth(8.8)

    ax1.errorbar(x, 2*meanRadius, stdRadius*2, marker='s', fillstyle='full')
    ax1.set_yscale('log')
    ax1.set_title('Diameter vs. Strahler order')
    ax2.errorbar(x, meanLength, stdLength, marker='s')
    ax2.set_yscale('log')
    ax2.set_title('Length vs. Strahler order')
    ax3.errorbar(x, meanVolume, stdVolume, marker='s')
    ax3.set_yscale('log')
    ax3.set_title('Volume vs. Strahler order')
    ax4.errorbar(x, meanAspectRatio, stdAspectRatio, marker='s')
    ax4.set_yscale('log')
    ax4.set_title('Aspect ratio vs. Strahler order')

    if outputImageFile:
        print('\nWriting outputs (image+data) in', outputImageFile[:-4] + '.png/.dat')
        data = np.vstack((x, meanRadius, stdRadius, meanLength, stdLength, meanVolume, stdVolume, meanAspectRatio, stdAspectRatio, orderDistribution))
        np.savetxt(outputImageFile[:-4] + '.dat', data)
        plt.savefig(outputImageFile)
    plt.show()
    
    return

def PlotLineAndBars(x, yline, ybars, labels=[r'Radius ($\mu m$)', 'Number of vessels per stream order'], fontsize=14, logplot=False):

    plt.style.use('seaborn')
    width = 345

    tex_fonts = {
            # Use LaTeX to write all text
            "text.usetex": True,
            "font.family": "serif",
            # Use 10pt font in plots, to match 10pt font in document
            "axes.labelsize": fontsize,
            "font.size": fontsize,
            # Make the legend/label fonts a little smaller
            "legend.fontsize": fontsize-2,
            "xtick.labelsize": fontsize-2,
            "ytick.labelsize": fontsize-2
        }

    plt.rcParams.update(tex_fonts)
    
    plt.figure()
    plt.rcParams['font.size'] = fontsize
    plt.xlabel('Horton-Srahler stream order')
    plt.ylim(yline[0,:].min()*0.75, yline[0,:].max()*1.25)
    plt.ylabel(labels[0])
    if logplot:
        plt.yscale('log')
    plt.errorbar(x, yline[0,:], yline[1,:], capsize=4, color='black')

    ax2 = plt.twinx()
    ax2.set_ylim(ybars.min()*0.75, ybars.max()*1.25)
    ax2.set_ylabel(labels[1])
    ax2.bar(x, ybars, alpha=0.3, color='gray')
    plt.show()
    return
