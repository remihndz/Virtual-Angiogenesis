#!/usr/bin/python
import sys

import numpy as np
from matplotlib.pyplot import imread
import matplotlib.pyplot as plt
from utils import FractalDimensionForParaviewScreenShot, FractalDimensionMap
 

threshold = 0

D = []
print('Found', len(sys.argv[1:]), 'files')
for imgFile in sys.argv[1:]:
    print('Reading file', imgFile)
    img = imread(imgFile)
    D.append(FractalDimensionForParaviewScreenShot(img, threshold=threshold))

D = np.array(D)
print('Mean and standard deviation:', D.mean(), D.std())
np.savetxt('FractalDimensionOfGeneratedTrees.dat', D)

# LocalFractalDimensionMap(img, w=1)
# FractalDimensionMap(img, w=2)
# FractalDimensionMap(img, w=3)

