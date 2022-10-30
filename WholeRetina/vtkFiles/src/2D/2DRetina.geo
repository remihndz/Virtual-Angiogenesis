//+
SetFactory("OpenCASCADE");
Disk(1) = {0, 0, 0, 1.4, 1.4}; // In cm
Characteristic Length{ PointsOf{ Surface{1}; } } = 0.01;