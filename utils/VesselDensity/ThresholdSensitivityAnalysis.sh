#!/bin/bash

outputfile="SensitivityToThreshold.dat"
rm $outputfile
echo Filename VD Threshold  >> outputfile

for filename in /mnt/c/Users/rhernand/Desktop/OCTADataSets/ROSE-dataset/data/ROSE-1/SVC/test/img/*.tif;
do    
    threshold="0"
    maxThreshold="255"
    while [ $threshold -lt $maxThreshold ]
    do
	python3 VesselDensity.py $filename $threshold >> $outputfile
	threshold=$((threshold+5))
    done
done

for filename in /mnt/c/Users/rhernand/Desktop/OCTADataSets/ROSE-dataset/data/ROSE-1/SVC/train/img/*.tif;
do    
    threshold="0"
    maxThreshold="255"
    while [ $threshold -lt $maxThreshold ]
    do
	python3 VesselDensity.py $filename $threshold >> $outputfile
	threshold=$((threshold+5))
    done
done

# Test segmented images as well
outputfile="SensitivityToThreshold_segmentedImages.dat"
for filename in /mnt/c/Users/rhernand/Desktop/OCTADataSets/ROSE-dataset/data/ROSE-1/SVC/test/gt/*.tif;
do    
    threshold="0"
    maxThreshold="255"
    while [ $threshold -lt $maxThreshold ]
    do
	python3 VesselDensity.py $filename $threshold >> $outputfile
	threshold=$((threshold+5))
    done
done

for filename in /mnt/c/Users/rhernand/Desktop/OCTADataSets/ROSE-dataset/data/ROSE-1/SVC/train/gt/*.tif;
do    
    threshold="0"
    maxThreshold="255"
    while [ $threshold -lt $maxThreshold ]
    do
	python3 VesselDensity.py $filename $threshold >> $outputfile
	threshold=$((threshold+5))
    done
done
