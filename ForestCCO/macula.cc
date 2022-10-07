#include <cmath>
#include <iostream>
#include <string>

#include "CompetingOptimizedArterialTrees.h"
#include "ForestCcoInvasion.h"
#include "CubeFunction.h"
#include "CircleFunction.h"
#include "DiscFunction.h"
#include "DomainFile.h"
#include "Tree.h"
#include "TreeFile.h"
#include "SuperficialVascularPlexus.h"

using namespace std;

int main(int argc, char *argv[]) {
  /* Declare the variables: */
  int numberOfTrees = 2,
    numberOfConnections = 20,	// Useless?
    maximumNumberOfAttempts = 10,
    numberOfTerminals = 200;

  // Domain parameters
  std::string vtkFileName = "Geometry/SVP";
  double AreaFAZ = 0.3e-6;	// The area of the FAZ, in m^2, from literature
  double innerRadius = sqrt(AreaFAZ/M_PI), // The radius of the FAZ (a circle), in m
    outerRadius = 3e-3,		// The radius of the FOV, in m
    inletFlow = 3.611e-6,	// 3.61 muL/min = 3.1e-7 L/min
    outletPressure = 2.133e3,   // 2133.16 N/m^2 (Pa) = 16 mm Hg
    outputUnit = 1e6;	        // 1e6 for microns

  std::cout << "Generating vessels in an annulus with inner radius "
	    << innerRadius
	    << "m and outer radius "
	    << outerRadius
	    << "m. Output unit is the micron." << std::endl;
  
  // Generate a vtk file with 6 random seeds and 50000 random points
  SuperficialVascularPlexus(vtkFileName, innerRadius, outerRadius,
			    numberOfTrees, 100000);

  
  // Perfusion parameters
  double targetPerfArt = 1.0/numberOfTrees,
    targetPerfVein     = 1.0/numberOfTrees;

  // CCO parameters
  double targetPerfusionFlow[numberOfTrees],
    stageCoefficient = 0.1,
    radiusExpoent = 2.75,
    lengthExpoent = 1.0;

  for (int i=0; i<numberOfTrees; i++)
    targetPerfusionFlow[i] = targetPerfArt;
  
  Domain *domain;
  DomainFile *domainFile;
  TreeModel **trees;
  TreeFile *treeFile;
 
  /* Instantiate the DomainFile and TreeModel: */
  
  domainFile = new DomainFile(vtkFileName+".vtk",			      
			      new DiscFunction(2, innerRadius, outerRadius)
  );


  trees = new TreeModel *[numberOfTrees];
  for (int i = 0; i < numberOfTrees; i++)
    {
      trees[i] = new Tree(
			  domainFile->seed(i),
			  numberOfTerminals,
			  domainFile->dimension()
			  );
   
      /* Set the perfusion volume and the terminal pressure: */
      /* The total perfusion flow is 0.00000833 m^3/s = 500 ml/min */
      // trees[i]->setPerfusionFlow(targetPerfusionFlow[0] * 8.33e-6);

      if (i<floor(numberOfTrees/2))
	trees[i]->setPerfusionFlow(targetPerfArt * inletFlow);
      else
	trees[i]->setPerfusionFlow(targetPerfVein * inletFlow);
      
      trees[i]->setTerminalPressure(outletPressure);

    }
  /* Instantiate the ForestCcoInvasion object: */
  CompetingOptimizedArterialTrees *coat = new CompetingOptimizedArterialTrees(
    domainFile,
    trees,
    numberOfTrees,
    numberOfTerminals,
    stageCoefficient,
    targetPerfusionFlow,
    radiusExpoent,
    lengthExpoent,
    0.0,
    0.0,
    M_PI/4.0
  );
  // ForestCcoInvasion *coat = new ForestCcoInvasion(domainFile,
  // 			       trees,
  // 			       numberOfTrees,
  // 			       numberOfTerminals,
  // 			       0.1, // Invasion coefficient
  // 			       targetPerfusionFlow,
  // 			       radiusExpoent,
  // 			       lengthExpoent);
  
  std::cout << "Initialization successful." << std::endl;

  /* Grow the tree: */
  coat->grow();
  
  /* Save the tree file: */
  for (int i = 0; i < numberOfTrees; i++)
    {
      treeFile = new TreeFile(coat->tree(i));

      /* 
	 Set the unit to microns (1m = 1e6micron)
      */
      treeFile->tree()->setLengthUnit(1e6);
      treeFile->tree()->setRadiusUnit(1e6);

      treeFile->save("svp_"+std::to_string(i)+".vtk");

    }
  
  return 0;
}
