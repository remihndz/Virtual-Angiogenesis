/**
 * @file CapillaryPlexus.h
 * @author Remi Hernandez (remi.hernandez@liverpool.ac.uk)
 * @date 05-10-2022
 * @version 1.0
 */

#ifndef CAPILLARY_PLEXUS_H
#define CAPILLARY_PLEXUS_H

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>

#include <random>
#include <iostream>
#include <fstream>
#include "Point.h"


 /**
  * @brief Construct a region to be populated by vessels and save it as a .vtk file.
  *
  * @param VTKFilename The name of the vtk file to write the output in.
  * @param InnerRadius The radius of the FAZ.
  * @param OuterRadius The radius of the FOV.
  * @param Thickness   Thickness of the complex.
  * @param Seeds       An array of seeds.
  * @param numberOfPoints The number of points to be generated.
  */
void CapillaryPlexus(std::string VTKFilename,
		     double InnerRadius,
		     double OuterRadius,
		     double Thickness,
		     std::vector<Point> Seeds,
		     int numberOfPoints);

 /**
  * @brief Construct a region to be populated by vessels and save it as a .vtk file.
  *
  * @param VTKFilename The name of the vtk file to write the output in.
  * @param InnerRadius The radius of the FAZ.
  * @param OuterRadius The radius of the FOV.
  * @param Thickness   Thickness of the complex.
  * @param Seeds       An array of seeds.
  */
void CapillaryPlexus(std::string VTKFilename,
		     double InnerRadius,
		     double OuterRadius,
		     double Thickness,
		     std::vector<Point> Seeds);




#endif	 // CAPILLARY_PLEXUS_H