import sys
from pathlib import Path
from utils import *
from numpy import loadtxt

# How many order to compute (stacks all vessels above that into the same order)
orderMax = 10



# Parameters of the plots
fontsize = 20
markersize = 10

# plt.style.use('classic')
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

treesData = []

bifRatios = [[],[]]
maxR, minR = [], []
for ccoFile in sys.argv[1:]:
    print("Reading file", Path(ccoFile).stem + '.cco')
    data, connectivity = ReadTree(ccoFile)

    orderedData = BifurcationOrder(data, connectivity, orderMax=orderMax)
    treesData.extend(orderedData)
    a,b = BifurcationDiameterRatio(data, connectivity, plot=False)

    bifRatios[0].extend(a)
    bifRatios[1].extend(b)

    rad = np.array(orderedData)[:,1]*1e3
    minR.append(rad.min())
    maxR.append(rad.max())

minR = np.array(minR)
maxR = np.array(maxR)
maxOrder = np.array(treesData)[:,-1].astype(int).max()
data = np.array(treesData)

dataRadius = np.zeros((maxOrder, len(treesData)))
dataLength = np.zeros((maxOrder, len(treesData)))
dataAspectRatio   = np.zeros((maxOrder, len(treesData)))
orderDistribution = np.zeros((maxOrder,))

# Load the data
for Id, radius, length, flow, stage, order in treesData:
    r,l = radius*1e3, length*1e3 # Convert to microns
    i,j = int(order)-1, int(Id)
    dataRadius[i, j] = r
    dataLength[i, j] = l
    dataAspectRatio[i, j] = l/r
    orderDistribution[i] += 1

orderDistribution /= len(sys.argv[1:])
dataRadius = dataRadius.transpose()

dataLength = dataLength.transpose()
dataAspectRatio = dataAspectRatio.transpose()
# Compute descriptive statistics for each order

meanRadius = np.mean(dataRadius, axis = 0, where=dataRadius>0)
meanLength = np.mean(dataLength, axis = 0, where=dataRadius>0)
meanAspectRatio = np.mean(dataAspectRatio, axis = 0, where=dataRadius>0)
minD, maxD = dataRadius[dataRadius>0].min()*2.0,  dataRadius.max()*2.0

# Compute std above and below mean
stdRadius = np.zeros((2, maxOrder))
stdLength = np.zeros_like(stdRadius)
stdAspectRatio = np.zeros_like(stdRadius)

countR, countL, countAsp =  np.ones((2, maxOrder)),  np.ones((2, maxOrder)),  np.ones((2, maxOrder))
for i in range(dataRadius.shape[0]):
    for j in range(dataRadius.shape[1]):
        meanR, meanL, meanAsp = meanRadius[j], meanLength[j], meanAspectRatio[j]
        if dataRadius[i,j] > 0:
            r, l, asp = dataRadius[i,j], dataLength[i,j], dataAspectRatio[i,j]
            kR, kL, kAsp = int(r>meanR), int(l<meanL), int(asp<meanAsp)
            stdRadius[kR,j] += (r-meanR)**2
            stdLength[kL,j] += (l-meanL)**2
            stdAspectRatio[kAsp,j] += (asp-meanAsp)**2
            countR[kR,j] += 1
            countL[kL,j] += 1
            countAsp[kAsp,j] += 1

stdRadius, stdLength, stdAspectRatio = np.divide(stdRadius, countR), np.divide(stdLength, countL), np.divide(stdAspectRatio, countAsp)
stdRadius, stdLength, stdAspectRatio = stdRadius**.5, stdLength**.5, stdAspectRatio**.5

cutoff = 9
meanRadiusBis, stdRadiusBis, orderDistributionBis  = np.zeros((cutoff+1,)), np.zeros((2,cutoff+1)), np.zeros((cutoff+1,))
meanRadiusBis[:-1] = meanRadius[:cutoff]
meanRadiusBis[-1] = np.mean(meanRadius[cutoff:])
stdRadiusBis[:,:-1] = stdRadius[:,:cutoff]
stdRadiusBis[:,-1] = np.nan_to_num(np.mean(stdRadius[cutoff:]**2)**.5)
orderDistributionBis[:-1] = orderDistribution[:cutoff]
orderDistributionBis[-1] = np.sum(orderDistribution[cutoff:])

meanRadius, stdRadius, orderDistribution = meanRadiusBis, stdRadiusBis, np.flip(orderDistributionBis)

# Load data from papers

x = np.arange(meanRadius.shape[0])
Takahashi = np.loadtxt('img/Takahashi.dat')[5:,:] # By column: order, diameter, length

# Uncomment for small to large classification
# Takahashi[:,0] = np.flip(Takahashi[:,0])-Takahashi[:,0].min()+1
An = np.loadtxt('img/An2020MeanDiameter.dat')   # By column: order, mean diameter

YuLength = np.loadtxt('img/Yu2010Length.dat')   # By column: order, mean length, std
YuDiameter = np.loadtxt('img/Yu2010Diameter.dat') # By column: order, mean diameter (width), std
KornfieldDiameter = np.loadtxt('img/Kornfield2014.dat')

plt.figure(1, figsize=[8.8, 7.2])
plt.xlabel('Number of bifurcations')
plt.ylabel(r'Diameter ($\mu m$)')

plt.errorbar(x, 2*meanRadius, 2*stdRadius, capsize=4, color='black',
             label=f'This work (diameter range: {(2*minR).mean():1.1f} to {(2*maxR).mean():1.1f})',
             marker='s', markersize=markersize)

plt.plot(np.flip(Takahashi[:,1]), label="Takahashi's ideal network", linestyle='-.', color='black')

plt.plot(x[-5:], np.flip(An[:,1]), label='An 2020, mean value', linestyle='--', marker='^', color='black', markersize=markersize)

plt.plot(x[-4:], np.flip(KornfieldDiameter[:,1]), label='Kornfield 2014, rat retina', linestyle='dotted', marker='v', color='black', markersize=markersize)

plt.legend()
ax2 = plt.twinx()
ax2.set_ylabel('Average number of vessels per group')

ax2.bar(x, np.flip(orderDistribution), alpha=0.3, color='gray')

print('Min/max diameter for my networks:', minD, maxD)
print('Min/max for Takahashi:', Takahashi[:,1].min(), ' / ', Takahashi[:,1].max())

ticks_labels = [str(xi) for xi in x[:-1]]
ticks_labels.append(r'$>$' + str(x[-2]))
plt.xticks(x, ticks_labels)

tableToSave = np.zeros((x.size, 5))
tableToSave[:,0] = x
tableToSave[:,1] = meanRadius
tableToSave[:,2] = stdRadius[0,:]
tableToSave[:,3] = stdRadius[1,:]
tableToSave[:,4] = np.flip(orderDistribution)
np.savetxt('RadiusPerBifurcationOrder.dat', tableToSave, header="BifOrder meanRadius AboveStd BelowStd Freq", comments='')
plt.savefig('DiameterVsBranchingOrder.png', dpi=500, bbox_inches='tight')


plt.figure(2, figsize=[7.2, 5.6])
bifRatios.sort()
bifRatios = np.array(bifRatios)

dsdl, dldp = np.array(bifRatios[0]), np.array(bifRatios[1])
tableToSave = np.zeros((dsdl.shape[0], 2))
tableToSave[:,0] = dsdl
tableToSave[:,1] = dldp
np.savetxt('DiameterRatiosAtBifurcations.dat', tableToSave, header="SmallerBranchOverLarger LargerBranchOverParent", comments='')

aggregatedData = []
x = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

for i,xi in enumerate(x[1:]):
    xim = x[i]
    mask = (dsdl>xim) & (dsdl<xi)
    aggregatedData.append(dldp[mask])
    
plt.boxplot(aggregatedData, positions=x[1:])

m = 2.85
x = np.linspace(0.2,1,500)
y = np.linspace(0.2,1,500)
X,Y = np.meshgrid(x,y)

f,g = Y**(-m), 1+X**m
z = f-g
cs1 = plt.contour(X,Y, z, [0], linestyles='--', linewidths=3)
cs1.collections[0].set_label(r"Murray's law with $\gamma=2.85$")

plt.xlabel("ds/dl")
plt.ylabel("dl/dp")
plt.xlim(0,1.2)
plt.ylim(0,1.2)
plt.legend()

plt.savefig('BifurcationRatio.png', dpi=300, bbox_inches='tight')

plt.show()

# Save all the data (large file)
f = open("TreesData.dat", "w")
f.write("Diameter Length Flow Stage Order\n")
for Id, radius, length, flow, stage, order in treesData:
        d,l = 2*radius*1e3, length*1e3 # Convert to microns
        order = min(order, orderMax)-1
        f.write(f"{d} {l} {flow} {stage} {order}\n")

f.close()


