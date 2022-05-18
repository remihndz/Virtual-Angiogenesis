import sys
from PIL import Image, ImageOps
import numpy as np
import matplotlib.pyplot as plt
from skimage.morphology import thin
from skimage.filters import roberts

file_names = sys.argv[1:]

Metrics = []

for file_name in file_names:

    img = ImageOps.grayscale(Image.open(file_name))

    VesselMap    = np.array(img.convert('1'))
    LengthMap    = thin(VesselMap)
    PerimeterMap = roberts(VesselMap)

    if len(file_names)==1:
        fig, ax = plt.subplots(2,2)
        ax1, ax2, ax3, ax4 = ax.ravel()
        ax1.imshow(img, cmap=plt.cm.gray)
        ax1.set_title('Original')

        ax2.imshow(VesselMap, cmap=plt.cm.gray)
        ax2.set_title('Vessel map')

        ax3.imshow(LengthMap, cmap=plt.cm.gray)
        ax3.set_title('Length map')

        ax4.imshow(PerimeterMap, cmap=plt.cm.gray)
        ax4.set_title('Perimeter map')

        for a in ax.ravel():
            a.axis('off')
        fig.tight_layout()
        plt.show()


    A,S,P,X = np.sum(VesselMap), np.sum(LengthMap), np.sum(PerimeterMap), VesselMap.size
    VAD = A/X
    VSD = S/X
    VCI = P**2/(4*np.pi*A)
    VDI = A/S
    VPI = P/X

    Metrics.append([VDI, VAD, VSD, VPI, VCI])

Metrics = np.array(Metrics)
np.savetxt('MetricsFromPNG.dat', Metrics)

    
ChuMetrics = [1, 0.518, 0.146, 0.349, 1]
Metrics[:,0] = Metrics[:,0]/24.421
Metrics[:,-1] = Metrics[:,-1]/17838.540
means = Metrics.mean(axis=0)
std = Metrics.std(axis=0)

labels = ['VDI (normalized)', 'VAD', 'VSD', 'VPI', 'VCI (normalized)']
bar_width=0.35
index = np.arange(5)

fig, ax = plt.subplots()
me = ax.bar(index, means, bar_width, yerr=std, label='From generated images')
chu = ax.bar(index+bar_width, ChuMetrics, bar_width, label='From OCTA')

ax.set_xlabel("Vascular Indices")
ax.set_xticks(index+bar_width/2)
ax.set_xticklabels(labels)
ax.legend()
plt.show()
