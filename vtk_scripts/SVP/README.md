# Creating empty geometries to be perfused

## Description
This folder contains scripts that create empty geometries to be used by the VItA codes. The geometries are generated as using the VTK library, which is compiled by default alongside the VItA library (https://github.com/GonzaloMaso/VItA).
If the VTK library has not been installed (i.e., not added to the default search path for shared libraries) on your system, make sure the path to the VTK library is added to your system's path to avoid runtime errors.

## Compilation
If VTK is installed:
```
cmake .
make
```
Alternatively, if VTK has only been compiled:
```
cmake -DVTK_DIR:PATH=/home/me/vtk_build .
make
```
