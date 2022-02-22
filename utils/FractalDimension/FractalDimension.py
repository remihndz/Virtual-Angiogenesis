#!/usr/bin/python
import sys

import numpy as np
from matplotlib.pyplot import imread
import matplotlib.pyplot as plt

def fractal_dimension(Z, threshold=10):
    # Only for 2d image
    assert(len(Z.shape) == 2)
    
    # From https://github.com/rougier/numpy-100 (#87)
    def boxcount(Z, k):
        S = np.add.reduceat(
            np.add.reduceat(Z, np.arange(0, Z.shape[0], k), axis=0),
            np.arange(0, Z.shape[1], k), axis=1)
        
        # We count non-empty (0) and non-full boxes (k*k)
        return len(np.where((S > 0) & (S < k*k))[0])
    
    # Transform Z into a binary array
    Z = (Z < threshold)
    plt.figure(1)
    plt.imshow(Z, cmap=plt.cm.gray_r)
    plt.title(f'Binarized image with threshold {threshold:1.2f}')
    plt.legend()
    
    # Minimal dimension of image
    p = min(Z.shape)

    # Greatest power of 2 less than or equal to p
    n = 2**np.floor(np.log(p)/np.log(2))
    
    # Extract the exponent
    n = int(np.log(n)/np.log(2))
    
    # Build successive box sizes (from 2**n down to 2**1)
    sizes = 2**np.arange(n, 1, -1)
    
    # Actual box counting with decreasing size
    counts = []
    for size in sizes:
        counts.append(boxcount(Z, size))
        
    # Fit the successive log(sizes) with log (counts)
    coeffs = np.polyfit(np.log(sizes), np.log(counts), 1)
    D = coeffs[0]
    
    print("Minkowski–Bouligand dimension (computed): ", -D)

    plt.figure(2)
    plt.plot(np.log(sizes), np.log(counts), '-s', label='Box count')
    plt.plot(np.log(sizes), np.log(counts[-1]) + D * (np.log(sizes)-np.log(sizes[-1])), label='Linear fit')
    plt.title(f'Minkowski-Bouligand dimension: {-D:1.2f}')
    plt.legend()
    plt.show()
    return


I = imread(str(sys.argv[1]))
if len(sys.argv)==3:
    threshold = float(sys.argv[2])
else:
    threshold = 50

fractal_dimension(I, threshold)
