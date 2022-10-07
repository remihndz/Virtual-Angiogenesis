import numpy as np
from numpy.random import random

n = 50000                          # Number of points in the domain (=possible distal ends for new vessels)
numberOfSeeds = 6                  # number of trees in the vasculature 
FOV       = 0.03                  # in m, the side of the square 
AreaFAZ   = 3.2e-6                 # in m^2, from https://doi.org/10.1371/journal.pone.0188572
RadiusFAZ = np.sqrt(AreaFAZ/np.pi) # in m

# Generate the points of the domain
R = FOV/np.sqrt(2)                          # The radius of the smallest circle enscribing the FOV: the boundary of the vascular domain
theta = random(n)*2.0*np.pi                 # Angle from the centre in [0, 2pi]
r     = RadiusFAZ + random(n)*(R-RadiusFAZ) # Distance from the centre
x,y = r*np.cos(theta), r*np.sin(theta)      # Position of the random points in Cartesian coordinates

seeds = random(numberOfSeeds)*2.0*np.pi           # Position on the outer circle
xSeeds, ySeeds = R*np.cos(seeds), R*np.sin(seeds) # Position of the seeds in Cartesian coordinates

print(f"Creating a {FOV}m wide slice of tissue, embedded in a circle of radius {R}m. The FAZ has radius {RadiusFAZ}m.")

with open('macula2D.vtk','w') as f:
    f.write("# vtk DataFile Version 3.0\n")
    f.write("Domain file for CCOLab 1.0\n")
    f.write("ASCII\n")
    f.write("FIELD macula 4\n")
    f.write("dimension 1 1 int\n")
    f.write("2\n\n")

    f.write("area 1 1 double\n")
    f.write(f"{R*R}\n\n")   # Area of the domain

    # Seed points
    f.write(f"seeds 2 {numberOfSeeds} double\n")
    for xi,yi in zip(xSeeds, ySeeds):
        f.write(f"{xi} {yi}\n")
    f.write("\n")

    # Random points within the domain
    f.write(f"points 2 {n} double\n")
    for xi, yi in zip(x,y):
        f.write(f"{xi} {yi}\n")
