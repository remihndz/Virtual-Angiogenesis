// STD Libs
#include<cmath>
#include<string>
#include<ctime>
#include<cstdio>
#include<cstdlib>
#include<vector>
#include<fstream>


// VItA Libs
#include<structures/tree/AbstractCostEstimator.h>
#include<structures/domain/AbstractDomain.h>
#include<structures/domain/StagedDomain.h>
#include<core/GeneratorData.h>
#include<core/StagedFRROTreeGenerator.h>
#include<constrains/ConstantConstraintFunction.h>
#include<structures/tree/SingleVesselCCOOTree.h>
#include<structures/tree/AbstractObjectCCOTree.h>
#include<structures/vascularElements/AbstractVascularElement.h>
#include<io/VTKObjectTreeNodalWriter.h>
#include<structures/domain/DomainNVR.h>
#include<structures/domain/SimpleDomain2D.h>
#include<structures/domain/NormalDistributionGenerator.h>
#include<structures/vascularElements/SingleVessel.h>
#include<structures/tree/AdimSproutingVolumetricCostEstimator.h>
#include<structures/tree/SproutingVolumetricCostEstimator.h>

#include<stats/ObjectTreeStatsManager.h>
#include<stats/VesselObjectHandler.h>

using namespace std;

double IntercapillaryDistance(SingleVesselCCOOTree *tree)
{
  double ICD {0.0};
  

  return ICD;
}


int main(int argc, char *argv[])
{
  // Read the configuration file passed as command line argument
  string testTreeCCOFile = argv[1];

  int nFail {200}, seed {2208}, nDraw {10000}, deltaNu {7};
  double lLimFr {0.9}, perfAreaFr {1.0}, closeNeighFr {4.0}, thetaMin {(3./18.)*M_PI};
  
  GeneratorData *genData = new GeneratorData(16000, nFail, lLimFr, perfAreaFr, closeNeighFr, 0.25, deltaNu, 0, false);
  AbstractConstraintFunction<double, int> *gam {new ConstantConstraintFunction<double, int>(2.85)};
  AbstractConstraintFunction<double, int> *delta {new ConstantConstraintFunction<double, int>(0.0)};
  AbstractConstraintFunction<double, int> *eta {new ConstantConstraintFunction<double, int>(3.6)};  
  SingleVesselCCOOTree *tree = new SingleVesselCCOOTree(testTreeCCOFile, genData, gam, delta, nu);

  double ICD = IntercapillaryDistance(tree);
  cout << "Tree's ICD is " << ICD << endl;
  
  return 0;
} 
