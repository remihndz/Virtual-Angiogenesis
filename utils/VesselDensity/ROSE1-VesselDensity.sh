#!/bin/bash

declare -i threshold=30
rm ROSE1_SVC_VD.dat ROSE1_SVC_DVC_VD.dat

shopt -s nullglob

echo Filename VesselDensity Threshold  >> ROSE1_SVC_VD.dat
for filename in /mnt/c/Users/rhernand/Desktop/OCTADataSets/ROSE-dataset/data/ROSE-1/SVC/test/img/*.tif; do
    python3 OCTAImageVesselDensity.py $filename $threshold >> ROSE1_SVC_VD.dat
done

for filename in /mnt/c/Users/rhernand/Desktop/OCTADataSets/ROSE-dataset/data/ROSE-1/SVC/train/img/*.tif; do
    python3 OCTAImageVesselDensity.py $filename $threshold >> ROSE1_SVC_VD.dat
done

# For comparison, compute vessel density on the combined SVC-DVC complex
echo Filename VesselDensity Threshold  >> ROSE1_SVC_DVC_VD.dat
for filename in /mnt/c/Users/rhernand/Desktop/OCTADataSets/ROSE-dataset/data/ROSE-1/SVC_DVC/test/img/*.png; do
    python3 OCTAImageVesselDensity.py $filename $threshold >> ROSE1_SVC_DVC_VD.dat
done

for filename in /mnt/c/Users/rhernand/Desktop/OCTADataSets/ROSE-dataset/data/ROSE-1/SVC_DVC/train/img/*.png; do
    python3 OCTAImageVesselDensity.py $filename $threshold >> ROSE1_SVC_DVC_VD.dat
done
shopt -u nullglob
