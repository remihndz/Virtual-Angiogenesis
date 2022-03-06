# Libs for image processing
from PIL import Image, ImageOps
from scipy.ndimage import distance_transform_edt
import cv2 as cv

# Other
import numpy as np
import matplotlib.pyplot as plt


def IntercapillaryDistanceMap(filename, plotResult=False):

    img = ImageOps.grayscale(Image.open(filename))
    # Crop the image to the edge of the FOV
    imgBox = img.getbbox()
    img = np.array(img.crop(imgBox))

    # Erode
    kernel = np.ones((2,2), np.uint8)
    erodedImg = cv.erode(img, kernel, iterations=1)

    # Compute distance map
    dHeight, dWidth = 3.0/img.shape[0], 3.0/img.shape[1] # Size of a pixel in mm
    D = distance_transform_edt(erodedImg==0, return_distances=True, sampling=[dHeight,dWidth])
    
    if plotResult:
        plt.subplot(121)
        plt.imshow(erodedImg)
        
        plt.subplot(122)
        plt.imshow(D, cmap='hot')
        plt.colorbar()
        plt.show()

    return D
