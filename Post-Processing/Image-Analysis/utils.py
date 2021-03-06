# Libs for image processing
from PIL import Image, ImageOps
from scipy.ndimage import distance_transform_edt
import cv2 as cv # Requires opencv-contrib-python for the thinning function


# Other
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm

def IntercapillaryDistanceMap(filename, plotResult=False):

    threshold = 0
    if filename[-4:] == '.tif':
        threshold = 50          # Threshold value that gives a VAD~0.5
        
    img = ImageOps.grayscale(Image.open(filename))
    img = img.point( lambda p: 255 if p > threshold else 0)
    img.convert('1')

    # Crop the image to the edges of the FOV
    imgBox = img.getbbox()
    img = np.array(img.crop(imgBox))

    # Erode
    # kernel = np.ones((2,2), np.uint8)
    # erodedImg = cv.erode(img, kernel, iterations=1)

    # Skeletonize (Zheng Suen algorithm)
    erodedImg = cv.ximgproc.thinning(img)

    # Compute distance map
    dHeight, dWidth = 3000/img.shape[0], 3000/img.shape[1] # Size of a pixel in microns
    D = distance_transform_edt(erodedImg==0, return_distances=True, sampling=[dHeight,dWidth])
    
    if plotResult:
        plt.subplot(121)
        plt.imshow(erodedImg, cmap='gray')
        
        plt.subplot(122)
        plt.imshow(D, cmap='viridis', norm=colors.PowerNorm(gamma=0.4, vmin=0.0, vmax=0.15, clip=True))
        plt.colorbar()
        plt.show()
        
    return D

def VesselDensityForOCTA(filename, plotResult=False, threshold=50):
        
    img = ImageOps.grayscale(Image.open(filename))
    img = img.point( lambda p: 1 if p > threshold else 0)
    img.convert('1')
    # # Crop the image to the edges of the FOV
    imgBox = img.getbbox()
    img = np.array(img.crop(imgBox))

    # # Erode
    kernel = np.ones((2,2), np.uint8)
    erodedImg = cv.erode(img, kernel, iterations=1)

    # Compute VAD

    if plotResult:
        plt.imshow(img, cmap='gray')
        plt.show()
    pix = np.array(img)         # Convert to an array
    VAD = pix.mean()

    return VAD




def FractalDimensionForParaviewScreenShot(img, threshold=0.0, plot=False):
    
    # Only for grayscale images
    if len(img.shape)>2:
        R, G, B = img[:,:,0], img[:,:,1], img[:,:,2]
        imgGray = 0.2989 * R + 0.5870 * G + 0.1140 * B
        img = imgGray
        

    # From https://github.com/rougier/numpy-100 (#87)
    def boxcount(img, k):
        S = np.add.reduceat(
            np.add.reduceat(img, np.arange(0, img.shape[0], k), axis=0),
            np.arange(0, img.shape[1], k), axis=1)
        
        # We count non-empty (0) and non-full boxes (k*k)
        return len(np.where((S > 0) & (S < k*k))[0])
    
    # Transform img into a binary array
    img = (img > threshold)
    
    # Minimal dimension of image
    p = min(img.shape)

    # Greatest power of 2 less than or equal to p
    n = 2**np.floor(np.log(p)/np.log(2))
    
    # Extract the exponent
    n = int(np.log(n)/np.log(2))
    
    # Build successive box sizes (from 2**n down to 2**1)
    sizes = 2**np.arange(n, 1, -1)
    
    # Actual box counting with decreasing size
    counts = []
    for size in sizes:
        counts.append(boxcount(img, size))
        
    # Fit the successive log(sizes) with log (counts)
    coeffs = np.polyfit(np.log(sizes), np.log(counts), 1)
    D = coeffs[0]
    
    if plot:
        print("Minkowski???Bouligand dimension (computed): ", -D)
            
        plt.figure(1)
        plt.imshow(img, cmap=plt.cm.gray_r)
        plt.title(f'Binarized image with threshold {threshold:1.2f}')
        plt.legend()

        plt.figure(2)
        plt.plot(np.log(sizes), np.log(counts), '-s', label='Box count')
        plt.plot(np.log(sizes), np.log(counts[-1]) + D * (np.log(sizes)-np.log(sizes[-1])), label='Linear fit')
        plt.title(f'Minkowski-Bouligand dimension: {-D:1.2f}')
        plt.legend()
        plt.show()

    return -D



def FractalDimension(img, threshold=15.0):

    # Only for grayscale images
    if len(img.shape)>2:
        R, G, B = img[:,:,0], img[:,:,1], img[:,:,2]
        imgGray = 0.2989 * R + 0.5870 * G + 0.1140 * B
        img = imgGray
        

    # From https://github.com/rougier/numpy-100 (#87)
    def boxcount(img, k):
        S = np.add.reduceat(
            np.add.reduceat(img, np.arange(0, img.shape[0], k), axis=0),
            np.arange(0, img.shape[1], k), axis=1)
        
        # We count non-empty (0) and non-full boxes (k*k)
        return len(np.where((S > 0) & (S < k*k))[0])
    
    # Transform img into a binary array
    img = (img < threshold)
    plt.figure(1)
    plt.imshow(img, cmap=plt.cm.gray_r)
    plt.title(f'Binarized image with threshold {threshold:1.2f}')
    plt.legend()
    
    # Minimal dimension of image
    p = min(img.shape)

    # Greatest power of 2 less than or equal to p
    n = 2**np.floor(np.log(p)/np.log(2))
    
    # Extract the exponent
    n = int(np.log(n)/np.log(2))
    
    # Build successive box sizes (from 2**n down to 2**1)
    sizes = 2**np.arange(n, 1, -1)
    
    # Actual box counting with decreasing size
    counts = []
    for size in sizes:
        counts.append(boxcount(img, size))
        
    # Fit the successive log(sizes) with log (counts)
    coeffs = np.polyfit(np.log(sizes), np.log(counts), 1)
    D = coeffs[0]
    
    print("Minkowski???Bouligand dimension (computed): ", -D)

    plt.figure(2)
    plt.plot(np.log(sizes), np.log(counts), '-s', label='Box count')
    plt.plot(np.log(sizes), np.log(counts[-1]) + D * (np.log(sizes)-np.log(sizes[-1])), label='Linear fit')
    plt.title(f'Minkowski-Bouligand dimension: {-D:1.2f}')
    plt.legend()
    plt.show()
    return -D


def FractalDimensionMap(img, w, threshold=0.0):
    
    # Only for grayscale images
    if len(img.shape)>2:
        R, G, B = img[:,:,0], img[:,:,1], img[:,:,2]
        imgGray = 0.2989 * R + 0.5870 * G + 0.1140 * B
        img = imgGray

    # Transform img into a binary array
    img = (img > threshold).astype(int) # Vessels in black, background in white
    # From https://github.com/rougier/numpy-100 (#87)
    def boxcount(img, k):
        S = np.add.reduceat(
            np.add.reduceat(img, np.arange(0, img.shape[0], k), axis=0),
            np.arange(0, img.shape[1], k), axis=1)
        
        # We count non-empty (0) and non-full boxes (k*k)
        return len(np.where((S > 0) & (S < k*k))[0])

    def LocalFractalDimension(localImg):
        if not np.any(localImg):
            return 0
        # Minimal dimension of image
        p = min(localImg.shape)
    
        # Greatest power of 2 less than or equal to p
        n = 2**np.floor(np.log(p)/np.log(2))

        # Extract the exponent
        n = int(np.log(n)/np.log(2))
        
        # Build successive box sizes (from 2**n down to 2**1)
        sizes = 2**np.arange(n, 0, -1)

        # Actual box counting with decreasing size
        counts = []
        for size in sizes:
            counts.append(boxcount(localImg, size))
        
        # Fit the successive log(sizes) with log (counts)
        coeffs = np.polyfit(np.log(sizes), np.log(counts), 1)
        D = coeffs[0]
        if np.isnan(D):
            return 0
        if D > 0:
            return 0

        return -D

    
    # Compute local map for each pixel
    FDMap = np.zeros_like(img)
    paddedImg = np.pad(img, ((w, w),), mode='constant', constant_values=0)
    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            k,l = i+w, j+w
            localImg = paddedImg[k-w:k+w+1, l-w:l+w+1]
            FDMap[i,j] = LocalFractalDimension(localImg)

    plt.subplot(121)
    plt.imshow(FDMap)
    plt.colorbar()
    plt.subplot(122)
    plt.imshow(img)
    plt.show()
    return FDMap
