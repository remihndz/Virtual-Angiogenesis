// Gmsh project created on Thu Oct 13 11:07:12 2022
SetFactory("OpenCASCADE");

Sphere(1) = {0, 0, 0, 1, -Pi/2, Pi/2-60*Pi/180, 2*Pi};
Characteristic Length{ PointsOf{ Volume{1}; } } = 0.008;

