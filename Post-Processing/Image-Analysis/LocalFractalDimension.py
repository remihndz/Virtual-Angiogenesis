#!/usr/bin/python
import sys

import numpy as np
from matplotlib.pyplot import imread
import matplotlib.pyplot as plt
from utils import FractalDimensionMap
import PIL.Image  


threshold = 0

D = []
for imgFile in sys.argv[1:]:
    img = PIL.Image.open(imgFile)
    imgBox = img.getbbox()
    img = img.crop(imgBox)
    img = np.array(img)
    FDMap = FractalDimensionMap(img, 3, threshold=threshold)
    for i in range(FDMap.shape[0]):
        for j in range(FDMap.shape[1]):
            if FDMap[i,j] < 0.3:
                FDMap[i,j] = 0.0
            elif FDMap[i,j] < 0.7:
                FDMap[i,j] = 0.7
            else:
                FDMap[i,j] = 1.0
    plt.imshow(FDMap)
    plt.colorbar()
    plt.show()
    D.append(FDMap.mean())


D = np.array(D)
print('Mean and standard deviation:', D.mean(), D.std())