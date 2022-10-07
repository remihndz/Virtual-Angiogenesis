/**
 * @file SuperficialVascularPlexus.cc
 * @author Remi Hernandez (remi.hernandez@liverpool.ac.uk)
 * @date 05-10-2022
 * @version 1.0
 */
#include "SuperficialVascularPlexus.h"

void SuperficialVascularPlexus(std::string VTKFilename,
			       double InnerRadius,
			       double OuterRadius,
			       int    numberOfSeeds,
			       int    numberOfPoints) {

  double theta,
    radius,
    x[numberOfPoints],
    y[numberOfPoints],
    xSeeds[numberOfSeeds],
    ySeeds[numberOfSeeds],
    thickness = InnerRadius;

  
  std::random_device randomDevice;
  std::mt19937 generator(randomDevice());
  std::uniform_real_distribution<double> thetaDistribution(0.0, 2.0*M_PI),
    radiusDistribution(InnerRadius, OuterRadius),
    thicknessDistribution(0.0, -thickness);

  // Generate the random numbers
  for (int i=0; i<numberOfPoints; i++) {

    theta  = thetaDistribution(generator);
    radius = radiusDistribution(generator);
    x[i]   = radius * cos(theta);
    y[i]   = radius * sin(theta);
  }

  for (int i=0; i<numberOfSeeds; i++) {

    theta = thetaDistribution(generator);
    xSeeds[i] = OuterRadius * cos(theta);
    ySeeds[i] = OuterRadius * sin(theta);
  }

  // Write the VTK file
  std::ofstream f;
  f.open(VTKFilename+".vtk");

  f << "# vtk DataFile Version 3.0" << std::endl;
  f << "Domain file for CCOLab 1.0" << std::endl;
  f << "ASCII" << std::endl;
  f << "FIELD macula 4" << std::endl;
  f << "dimension 1 1 int" << std::endl;
  f << "2" << std::endl << std::endl;

  f << "area 1 1 double" << std::endl;
  f << M_PI * OuterRadius << std::endl << std::endl; // Area of the domain

  // Seed points
  f << "seeds 3 " << numberOfSeeds << " double" << std::endl;
  for (int i = 0; i < numberOfSeeds; i++) 
    f << xSeeds[i] << " "
      << ySeeds[i] << " "
      << thicknessDistribution(generator) << std::endl;
  f << std::endl;
  
  
  // Random points within the domain
  f << "points 3 " << numberOfPoints << " double" << std::endl;
  for (int i = 0; i < numberOfPoints; i++)
    f << x[i] << " "
      << y[i] << " "
      << thicknessDistribution(generator) << std::endl;

  f.close();
}  
  
void SuperficialVascularPlexus(std::string VTKFilename,
			       double InnerRadius,
			       double OuterRadius) {
  SuperficialVascularPlexus(VTKFilename,
			    InnerRadius,
			    OuterRadius,
			    6,	// Number of seeds
			    50000); // Number of points
}

void SuperficialVascularPlexus(std::string VTKFilename,
			       double InnerRadius,
			       double OuterRadius,
			       int numberOfSeeds) {
  SuperficialVascularPlexus(VTKFilename,
			    InnerRadius,
			    OuterRadius,
			    numberOfSeeds,
			    50000); // Number of points
}
