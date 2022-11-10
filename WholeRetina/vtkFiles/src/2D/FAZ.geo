//+
SetFactory("OpenCASCADE");
Disk(1) = {0.5, 0, 0, 0.027, 0.025};
Characteristic Length{ PointsOf{ Surface{1}; } } = 0.01;