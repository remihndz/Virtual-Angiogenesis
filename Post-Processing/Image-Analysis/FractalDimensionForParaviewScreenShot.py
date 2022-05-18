#!/usr/bin/python
import sys

import numpy as np
from matplotlib.pyplot import imread
import matplotlib.pyplot as plt
from utils import FractalDimensionForParaviewScreenShot, FractalDimensionMap

 

threshold = 0

D = []
print('Found', len(sys.argv[1:]), 'files')

f = open('FD.dat', 'w')
f.write("File FD\n")

for imgFile in sys.argv[1:]:
    img = imread(imgFile)
    D.append(FractalDimensionForParaviewScreenShot(img, threshold=threshold))

    print('For', imgFile, 'FD=', D[-1])
    f.write(f"{imgFile} {D[-1]}\n")
    
D = np.array(D)
print('Mean and standard deviation:', D.mean(), D.std())