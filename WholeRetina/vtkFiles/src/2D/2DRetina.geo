//+
SetFactory("OpenCASCADE");
Disk(1) = {0, 0, 0, 1.15, 1.15}; // In cm
Characteristic Length{ PointsOf{ Surface{1}; } } = 0.01;