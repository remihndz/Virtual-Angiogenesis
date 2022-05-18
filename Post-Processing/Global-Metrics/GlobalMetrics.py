import sys
from utils import *
from math import pi

assert(len(sys.argv) > 1)
listOfFiles = sys.argv[1:]

domainAreaInCm = 0.09 - 0.04**2*pi

VAD, VPI, VDI, VCI, VSD, VBC = StatisticsMultipleTrees(listOfFiles, domainAreaInCm=domainAreaInCm)

print(f"Found {len(listOfFiles):d} files.")
print("Metric  |  Mean+-std")
print(f"VAD     |  {VAD.mean()} +- {VAD.std()}")
print(f"VPI     |  {VPI.mean()} +- {VPI.std()}")
print(f"VDI     |  {VDI.mean()} +- {VDI.std()}")
print(f"VCI     |  {VCI.mean()} +- {VCI.std()}")
print(f"VSD     |  {VSD.mean()} +- {VSD.std()}")
print(f"VBC     |  {VBC.mean()} +- {VBC.std()}")


f = open("GlobalMetrics.dat", "w")
f.write("File VAD VPI VDI VCI VSD VBC\n")

for file, vad, vpi, vdi, vci, vsd, vbc in zip(listOfFiles, VAD, VPI, VDI, VCI, VSD, VBC):
    f.write(f"{file} {vad} {vpi} {vdi} {vci} {vsd} {vbc}\n")

f.close()

f = open("GlobalMetrics_forPandas.dat", "w")
f.write("Patient Metric Value\n")
    
for file, vad, vpi, vdi, vci, vsd, vbc in zip(listOfFiles, VAD, VPI, VDI, VCI, VSD, VBC):
    patient=file[33]
    f.write(f"{patient} VAD {vad}\n{patient} VPI {vpi}\n{patient} VDI {vdi}\n{patient} VCI {vci}\n{patient} VSD {vsd}\n{patient} VBC {vbc}\n")

f.close()