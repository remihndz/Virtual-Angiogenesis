import numpy as np
import csv
import matplotlib.pyplot as plt

with open('SensitivityToThreshold.dat') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=' ')
    raw_data = []
    for row in csv_reader:
        raw_data.append(row)

count = 0
current_file = raw_data[0][0]
old_file = current_file
while current_file == old_file:
    count+=1
    current_file = raw_data[count][0]

nbExperiments = int(len(raw_data)/count)


vascularDensity = np.array([vd for misc,vd,thresh in raw_data]).reshape((nbExperiments, -1)).astype(float)
threshold = np.array([thresh for misc, vd, thresh in raw_data])[:count].astype(float)

# Compute mean and std
mean = np.mean(vascularDensity, axis=0)
std  = np.std(vascularDensity, axis=0)

plt.plot(threshold, mean, color='blue')
plt.fill_between(threshold, mean-std, mean+std, color='blue', alpha=0.5)




with open('SensitivityToThreshold_segmentedImages.dat') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=' ')
    raw_data = []
    for row in csv_reader:
        raw_data.append(row)

count = 0
current_file = raw_data[0][0]
old_file = current_file
while current_file == old_file:
    count+=1
    current_file = raw_data[count][0]

nbExperiments = int(len(raw_data)/count)


vascularDensity = np.array([vd for misc,vd,thresh in raw_data]).reshape((nbExperiments, -1)).astype(float)
threshold = np.array([thresh for misc, vd, thresh in raw_data])[:count].astype(float)

# Compute mean and std
mean = np.mean(vascularDensity, axis=0)
std  = np.std(vascularDensity, axis=0)

plt.plot(threshold, mean, color='yellow')
plt.fill_between(threshold, mean-std, mean+std, color='yellow', alpha=0.5)

plt.show()
