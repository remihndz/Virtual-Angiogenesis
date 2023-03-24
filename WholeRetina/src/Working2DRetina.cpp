// STD Libs
#include<cmath>
#include<string>
#include<ctime>
#include<cstdio>
#include<cstdlib>
#include<vector>
#include<fstream>
#include<omp.h>

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
#include<creators/ParallelepipedCreator.h>
#include <creators/CylinderCreator.h>

using namespace std;
string rootTreeFilename = "CRA.cco";
double q0 {15.0 * 0.001/60.0};	// Converts muL/min to cm^3/s
double r0 {0.0007};		// mm
point x0 {0.0, 0.0, -0.1};	// mm
double p0 {50 * 133.3224};  // mmHg*133.3224 = Pa, should be consistent with pressure unit 

string outputFilename = "../Results/2DRetina",
  hullVTKFilename = "../vtkFiles/Results/2D/Hull.vtp",
  FAZFilename = "../vtkFiles/Results/2D/FAZ.vtp";


void Vascularise(string outputFilename,
		 vector<long long int> nTerms,
		 vector<AbstractConstraintFunction<double, int>*> gammas,
		 vector<AbstractConstraintFunction<double, int>*> deltas,
		 vector<AbstractConstraintFunction<double, int>*> etas,
		 int nDraw, int seed, int nFail, double lLimFactor,
		 double perfusionAreaFactor, double closeNeighborhoodFactor, int DeltaNu,
		 double thetaMin)
{
  double lb[3] = {-0.3, -0.3, -0.1};
  double lu[3] = {0.3, 0.3, 0.1};
  ParallelepipedCreator *creator = new ParallelepipedCreator(lb, lu);
  creator->create(hullVTKFilename);
  delete creator;
  CylinderCreator *creator1 = new CylinderCreator({0.15, 0.0, -0.1}, 0.025, 0.2, 100);
  creator1->create(FAZFilename);
  delete creator1;

  cout << "DeltaNu = " << DeltaNu << endl;
  GeneratorData *generatorData = new GeneratorData(16000, nFail, lLimFactor, perfusionAreaFactor, 0.25, 1.0, DeltaNu, 0, false);
  DomainNVR *domain = new DomainNVR(hullVTKFilename, {FAZFilename}, generatorData);
  domain->setMinBifurcationAngle(thetaMin);
  domain->setIsConvexDomain(true);

  StagedDomain *stagedDomain = new StagedDomain();
  stagedDomain->addStage(nTerms[0], domain);

  StagedFRROTreeGenerator *treeGenerator = new StagedFRROTreeGenerator(stagedDomain, x0, r0, q0, nTerms[0], gammas, deltas, etas, 0., 1e-5);
  SingleVesselCCOOTree *tree = static_cast<SingleVesselCCOOTree *>(treeGenerator->getTree());
  tree->setIsInCm(false);

  tree = (SingleVesselCCOOTree *)treeGenerator->generate(10, ".");
  tree->save(outputFilename+".cco");
  VTKObjectTreeNodalWriter *treeWriter = new VTKObjectTreeNodalWriter();
  treeWriter->write(outputFilename + ".vtp", tree);
    
}

int main(int argc, char *argv[])
{
#pragma omp parallel
  {
    cout << omp_get_num_threads() << " threads running." << endl;
  }


  vector<long long int> nTerms = {100};
  vector<AbstractConstraintFunction<double, int>*> gammas = {new ConstantConstraintFunction<double, int>(3.0)},
    deltas = {new ConstantConstraintFunction<double, int>(0.0)},
    etas = {new ConstantConstraintFunction<double, int>(3.6)};
  double thetaMin {0}, // In radian?
    perfusionAreaFactor {.9},
    closeNeighborhoodFactor {1.0},
    lLimFactor {0.9};
  
  // Consecutive attempts to generate a point - nFail
  int nFail = 200;
  // Discretisation of the testing triangle for the bifurcation - Delta nu - Figure 1
  int DeltaNu = 7;
  // Buffer size for random point generation
  int nDraw {10000};
  // Random seed
  long long int seed {time(nullptr)};
  

  Vascularise(outputFilename,
	      nTerms,
	      gammas,
	      deltas,
	      etas,
	      nDraw, seed, nFail, lLimFactor,
	      perfusionAreaFactor, closeNeighborhoodFactor, DeltaNu,
	      thetaMin);

  for (auto p : gammas)
    delete p;
  gammas.clear();
  for (auto p : etas)
    delete p;
  etas.clear();
  for (auto p : deltas)
    delete p;
  deltas.clear();
  return 0;
} 
