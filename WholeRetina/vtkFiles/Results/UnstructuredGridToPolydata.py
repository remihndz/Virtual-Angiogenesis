# coding: utf-8
import sys
import vtkmodules.all as vtk

for filename in sys.argv[1:]:
    geofilter = vtk.vtkGeometryFilter()
    polydata = vtk.vtkPolyData()
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()
    geofilter.UnstructuredGridExecute(data, polydata)
    polydata
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileVersion(42)
    writer.SetFileName(filename[:-4] +'.vtp')
    writer.SetInputData(polydata)
    writer.Update()
    writer.Write()
