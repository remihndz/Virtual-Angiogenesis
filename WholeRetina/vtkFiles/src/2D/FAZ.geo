//+
SetFactory("OpenCASCADE");
Disk(1) = {0.2, 0, 0, 0.05, 0.04};
Characteristic Length{ PointsOf{ Surface{1}; } } = 0.01;