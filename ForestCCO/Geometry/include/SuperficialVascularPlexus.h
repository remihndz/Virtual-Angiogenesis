/**
 * @file SuperficialVascularPlexus.h
 * @author Remi Hernandez (remi.hernandez@liverpool.ac.uk)
 * @date 05-10-2022
 * @version 1.0
 */

#ifndef SUPERFICIAL_VASCULAR_PLEXUS_H
#define SUPERFICIAL_VASCULAR_PLEXUS_H

#define _USE_MATH_DEFINES
#include <cmath>

#include <random>
#include <iostream>
#include <fstream>


 /**
  * @brief Construct a region to be populated by vessels and save it as a .vtk file.
  *
  * @param VTKFilename The name of the vtk file to write the output in.
  * @param InnerRadius The radius of the FAZ
  * @param OuterRadius The radius of the FOV
  * @param numberOfSeeds The number of seeds to be generated
  * @param numberOfPoints The number of points to be generated
  */
void SuperficialVascularPlexus(std::string VTKFilename,
			       double InnerRadius,
			       double OuterRadius,
			       int numberOfSeeds,
			       int numberOfPoints);

 /**
  * @brief Construct a region to be populated by vessels and save it as a .vtk file.
  *
  * @param VTKFilename The name of the vtk file to write the output in.
  * @param InnerRadius The radius of the FAZ
  * @param OuterRadius The radius of the FOV
  */
void SuperficialVascularPlexus(std::string VTKFilename,
			       double InnerRadius,
			       double OuterRadius);

 /**
  * @brief Construct a region to be populated by vessels and save it as a .vtk file.
  *
  * @param VTKFilename The name of the vtk file to write the output in.
  * @param InnerRadius The radius of the FAZ
  * @param OuterRadius The radius of the FOV
  * @param numberOfSeeds The number of seeds to be generated
  */
void SuperficialVascularPlexus(std::string VTKFilename,
			       double InnerRadius,
			       double OuterRadius,
			       int numberOfSeeds);


#endif	 // SUPERFICIAL_VASCULAR_PLEXUS_H
