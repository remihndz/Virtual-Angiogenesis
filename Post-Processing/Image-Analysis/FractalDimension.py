#!/usr/bin/python
import sys

import numpy as np
from matplotlib.pyplot import imread
import matplotlib.pyplot as plt
from utils import FractalDimension
 

threshold = 30
makePlots = False

D = []
for imgFile in sys.argv[1:]:
    img = imread(imgFile)
    D.append(FractalDimension(img, threshold=threshold, plot=makePlots))

D = np.array(D)
print('Mean and standard deviation:', D.mean(), D.std())
np.savetxt('FractalDimensionOfOCTAs.dat', D)
