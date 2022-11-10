//+
SetFactory("OpenCASCADE");
Disk(1) = {0.5, 0, 0, 0.3, 0.3}; // In cm
Characteristic Length{ PointsOf{ Surface{1}; } } = 0.01;