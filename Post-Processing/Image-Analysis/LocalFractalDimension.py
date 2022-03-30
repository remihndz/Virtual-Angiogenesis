#!/usr/bin/python
import sys

import numpy as np
from matplotlib.pyplot import imread
import matplotlib.pyplot as plt
from utils import FractalDimensionMap
 

threshold = 0

D = []
for imgFile in sys.argv[1:]:
    img = imread(imgFile)
    D.append(FractalDimensionMap(img, 5, threshold=threshold).mean())

D = np.array(D)
print('Mean and standard deviation:', D.mean(), D.std())


