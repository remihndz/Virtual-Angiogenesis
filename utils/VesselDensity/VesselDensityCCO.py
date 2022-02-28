import sys
from utils import *

assert(len(sys.argv) > 1)
listOfFiles = sys.argv[1:]

VAD, VPI, VDI, VCI = StatisticsMultipleTrees(listOfFiles)
print("Metric  |  Mean+-std")
print(f"VAD     !  {VAD.mean()} += {VAD.std()}")
print(f"VPI     !  {VPI.mean()} += {VPI.std()}")
print(f"VDI     !  {VDI.mean()} += {VDI.std()}")
print(f"VCI     !  {VCI.mean()} += {VCI.std()}") 

