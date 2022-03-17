# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 14:15:11 2022

@author: rhernand
"""

import vtkmodules.all as vtk

# Read the source file.
reader = vtk.vtkXMLPolyDataReader()
reader.SetFileName(r"C:\Users\rhernand\Desktop\MyCodes\Virtual-Angiogenesis\Results50Sims\sim1.vtp")
reader.Update()  # Needed because of GetScalarRange
polyData = reader.GetOutput()


# Clipping
boundingBox = vtk.vtkBox()
boundingBox.SetBounds(-0.15, 0.15,
                      -0.15, 0.15,
                      -0.15, 0.15)

clipper = vtk.vtkClipPolyData()
clipper.SetInputData(polyData)
clipper.SetClipFunction(boundingBox)
clipper.SetInsideOut(True)
clipper.Update()

polyData = clipper.GetOutput()


whiteImage = vtk.vtkImageData()
bounds = polyData.GetBounds()
spacing = [0.05, 0.05, 0.001]
dim = [60, 60, 1]
whiteImage.SetDimensions(dim)

# origin = [bounds[0]+spacing[0]/2.0,
#           bounds[2]+spacing[1]/2.0,
#           bounds[4]+spacing[2]/2.0]
whiteImage.SetOrigin(-0.15,-0.15,0)
whiteImage.AllocateScalars(vtk.VTK_UNSIGNED_CHAR, 1)
inval = 255
outval = 0
count = whiteImage.GetNumberOfPoints()
for i in range(count):
    point = whiteImage.GetPointData()
    scalar = point.GetScalars()
    scalar.SetTuple1(i, inval)


pd2stenc = vtk.vtkPolyDataToImageStencil()
# pd2stenc.SetTolerance(0.5)
pd2stenc.SetInputData(polyData)
pd2stenc.SetInformationInput(whiteImage)
# pd2stenc.SetOutputOrigin(origin)
# pd2stenc.SetOutputSpacing(spacing)
# pd2stenc.SetOutputWholeExtent(whiteImage.GetExtent())
pd2stenc.Update()

imgstenc = vtk.vtkImageStencil()
imgstenc.SetInputData(whiteImage)
imgstenc.SetStencilConnection(pd2stenc.GetOutputPort())
imgstenc.ReverseStencilOff()
imgstenc.SetBackgroundValue( inval )
imgstenc.Update()

img = vtk.vtkImageData()
img = imgstenc.GetOutput()




