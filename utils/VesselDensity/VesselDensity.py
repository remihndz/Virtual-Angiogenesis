#!/usr/bin/python
import sys

import numpy as np
from matplotlib.pyplot import imread
import matplotlib.pyplot as plt

def VesselDensity(img, threshold=20.0, plot=False):

    # Only for grayscale images
    if len(img.shape)>2:
        R, G, B = img[:,:,0], img[:,:,1], img[:,:,2]
        imgGray = 0.2989 * R + 0.5870 * G + 0.1140 * B
        img = imgGray*255.0
    
    ''' Transform into binary image using the
        specified threshold value. '''
    img_bin = (img > threshold).astype(float)

    if plot:
        fig1 = plt.figure(1)
        plt.subplot(121)
        plt.imshow(img_bin, cmap=plt.cm.gray_r)
        plt.title(f'Binarized imaged. Threshold = {threshold:2.1f}.')
        
        plt.subplot(122)
        plt.imshow(img, cmap=plt.cm.gray_r)
        plt.title('Original image.')
        plt.show()

    TotalArea     = img.shape[0]*img.shape[1]
    VesselArea    = np.sum(np.sum(img_bin))
    VesselDensity = VesselArea/TotalArea * 100.00

    print(f"{sys.argv[1]} {VesselDensity:2.2f} {threshold:2.2f}")
    
    
img = imread(str(sys.argv[1]))

if len(sys.argv)==3:
    threshold = float(sys.argv[2])
else:
    threshold = 50.0

VesselDensity(img, threshold=threshold, plot=False)

