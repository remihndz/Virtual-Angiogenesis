import sys
from pathlib import Path
from utils import *
from numpy import loadtxt

# Parameters of the plot
fontsize = 16
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
bifRatios = []
maxR, minR = [], []
for ccoFile in sys.argv[1:]:
    print("Reading file", Path(ccoFile).stem + '.cco')
    data, connectivity = ReadTree(ccoFile)

    orderedData = BifurcationOrder(data, connectivity)
    treesData.extend(orderedData)
    a,b = BifurcationDiameterRatio(data, connectivity, plot=False)
    bifRatios.append([a,b])
    rad = np.array(data)[:,1]*1e3
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
np.savetxt('./Radii.dat', dataRadius)
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

countR, countL, countAsp =  np.zeros((2, maxOrder)),  np.zeros((2, maxOrder)),  np.zeros((2, maxOrder))
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

# Load data from papers
x = np.arange(1, maxOrder+1)
Takahashi = np.loadtxt('img/Takahashi.dat')[3:,:] # By column: order, diameter, length
Takahashi[:,0] = np.flip(Takahashi[:,0])-Takahashi[:,0].min()+1
An = np.loadtxt('img/An2020MeanDiameter.dat')   # By column: order, mean diameter
# YuLength = np.loadtxt('img/Yu2010Length.dat')   # By column: order, mean length, std
# YuDiameter = np.loadtxt('img/Yu2010Diameter.dat') # By column: order, mean diameter (width), std
KornfieldDiameter = np.loadtxt('img/Kornfield2014.dat')

plt.figure(1)
plt.xlabel('Horton-Strahler stream order')
plt.ylabel(r'Diameter ($\mu m$)')
plt.errorbar(x, 2*meanRadius, 2*stdRadius, capsize=4, color='black',
             label=f'This work ({(2*minR).mean():1.2f}+-{(2*minR).std():1.2f}, {(2*maxR).mean():1.2f}+-{(2*maxR).std():1.2f})',
             marker='s')

# plt.plot(Takahashi[:,0], Takahashi[:,1], label="Takahashi's ideal network", linestyle='-.', color='black')

# plt.plot(x[-int(An[:,0].max()):], np.flip(An[:,1]), label='An 2020, mean value', linestyle='--', marker='^', color='black')

# plt.plot(x[-int(KornfieldDiameter[:,0].max()):], np.flip(KornfieldDiameter[:,1]), label='Kornfield 2014, rat retina', linestyle='dotted', marker='v', color='black')


# plt.errorbar(YuDiameter[:,0]+6, YuDiameter[:,1], YuDiameter[:,2], label=r'Yu 2020, mean$\pm$std', linestyle='dotted', marker='v', color='black', capsize=4)

# plt.yscale('log')
# plt.grid(True, which='both', linestyle='dotted')
plt.legend()
ax2 = plt.twinx()
ax2.set_ylabel('Average number of vessels per stream order')
ax2.bar(x, orderDistribution, alpha=0.3, color='gray')

# plt.figure(2)
# plt.xlabel('Horton-Strahler stream order')
# plt.ylabel(r'Length ($\mu m$)')
# plt.errorbar(x, meanLength, stdLength, capsize=4, color='black')
# plt.plot(Takahashi[:,0], Takahashi[:,2], label="Takahashi's ideal network", linestyle='-.', color='black')

# # plt.errorbar(YuDiameter[:,0]+6, YuLength[:,1], YuLength[:,2], label=r'Yu 2020, mean$\pm$std', marker='v', linestyle='dotted', color='black', capsize=4)

# plt.legend()
# ax2 = plt.twinx()
# ax2.set_ylabel('Average number of vessels per stream order')
# ax2.bar(x, orderDistribution, alpha=0.3, color='gray')

# plt.figure(3)
# plt.xlabel('Horton-Strahler stream order')
# plt.ylabel(r'Aspect ratio')
# plt.errorbar(x, meanAspectRatio, stdAspectRatio, capsize=4, color='black', label=r'This work, mean$\pm$std')
# plt.legend()
# ax2 = plt.twinx()
# ax2.set_ylabel('Average number of vessels per stream order')
# ax2.bar(x, orderDistribution, alpha=0.3, color='gray')

print('Min/max diameter for my networks:', minD, maxD)
print('Min/max for Takahashi:', Takahashi[:,1].min(), ' / ', Takahashi[:,1].max())


plt.figure(2)
bifRatios.sort()
bifRatios = np.array(bifRatios)
dsdl, dldp = bifRatios[:,0], bifRatios[:,1]

aggregatedData = []
x = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
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
m = 3
f,g = Y**(-m), 1+X**m
cs2 = plt.contour(X,Y,(f-g), [0], linestyles='-.', linewidths=3)

cs1.collections[0].set_label('Theoretical, m=2.85')
cs2.collections[0].set_label('Theoretical, m=3')

plt.xlabel("ds/dl")
plt.ylabel("dl/dp")
plt.xlim(0,1.2)
plt.ylim(0,1.2)
plt.legend()


plt.show()

