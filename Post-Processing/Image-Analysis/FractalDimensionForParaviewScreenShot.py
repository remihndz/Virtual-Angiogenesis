#!/usr/bin/python
import sys

import numpy as np
from matplotlib.pyplot import imread
import matplotlib.pyplot as plt
from utils import FractalDimensionForParaviewScreenShot, LocalFractalDimensionMap
 

threshold = 0

D = []
for imgFile in sys.argv[1:]:
    img = imread(imgFile)
    D.append(FractalDimensionForParaviewScreenShot(img, threshold=threshold))

D = np.array(D)
print('Mean and standard deviation:', D.mean(), D.std())
np.savetxt('FractalDimensionOfGeneratedTrees.dat', D)

# LocalFractalDimensionMap(img, w=1)
LocalFractalDimensionMap(img, w=2)
LocalFractalDimensionMap(img, w=3)

