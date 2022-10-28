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
#include<structures/tree/VolumetricCostEstimator.h>
#include<structures/tree/AdimSproutingVolumetricCostEstimator.h>
#include<structures/tree/SproutingVolumetricCostEstimator.h>
#include<creators/ParallelepipedCreator.h>


using namespace std;
string rootTreeFilename = "CRA.cco";
double q0 {15.0 * 0.001/60.0};	// Converts muL/min to cm^3/s
double r0 {0.07};		// mm
point x0 {0.26, 0.0, 0.0};	// mm
double p0 {50 * 133.3224};  // mmHg*133.3224 = Pa, should be consistent with pressure unit 

string outputFilename = "../Results/2DRetina",
  hullVTKFilename = "../vtkFiles/Results/2D/Hull.vtk";
vector<string> NVRVTKFilenames({"../vtkFiles/Results/2D/FAZ.vtk"});


void Vascularise(string outputFilename,
		 string hullVTKFilename,
		 vector<string> NVRVTKFilenames,
		 vector<long long int> nTerms,
		 vector<AbstractConstraintFunction<double, int>*> gammas,
		 vector<AbstractConstraintFunction<double, int>*> deltas,
		 vector<AbstractConstraintFunction<double, int>*> etas,
		 int nDraw, int seed, int nFail, double lLimFactor,
		 double perfusionAreaFactor, double closeNeighborhoodFactor, int DeltaNu,
		 double thetaMin)
{
  // Here create GeneratorDatas and Domains (one of each for each stage)
  // TODO: change the input parameters to vectors to have different parameters for each stage
  vector<GeneratorData*> generators;
  vector<DomainNVR*> domains;
  StagedDomain *stagedDomain = new StagedDomain();
  // Does initializing pointers in a loop work? Or are the objects created deleted after the loop?
  // for (int i=0; i<nTerms.size(); i++)
  //   {
  //     SproutingVolumetricCostEstimator *FSprout = new SproutingVolumetricCostEstimator(50, 0.5, 1e+4);
  //     AbstractCostEstimator *costEstimator = FSprout;

  //     generators.push_back(new GeneratorData(16000, // Levels for tree scaling for each new segment test.
  // 					     nFail, // Number of trials before diminish dlim.
  // 					     lLimFactor, // Factor by which the Dlim constraint diminish after N failed trials.
  // 					     perfusionAreaFactor, // Factor that scales the perfusion area by which Dlim is computed.
  // 					     closeNeighborhoodFactor, // Factor that increase the neighborhood to search nearest neighbors.
  // 					     0.1, // Factor to scale the dLim to the middle point of the new vessel to avoid close neighbors.
  // 					     DeltaNu, // Number of bifurcation sites tested in the optimization process is given by nBifurcationTest * ( nBifurcationTest - 1 ). (default 8)
  // 					     1,	   // Functionality of the vessel generated, important for Object trees.
  // 					     false, // Indicates if dLimCorrectionFactor must be resetted to 1 when the stage begins.
  // 					     costEstimator)); // Cost estimator for the given stage 

  //     domains.push_back(new DomainNVR(hullVTKFilename, NVRVTKFilenames,
  // 				      nDraw, seed, generators.back()));
  //     domains.back()->setIsConvexDomain(true);
  //     domains.back()->setMinBifurcationAngle(thetaMin);

  //     stagedDomain->addStage(nTerms[i], domains.back());
  //   }

  SproutingVolumetricCostEstimator *FSprout = new SproutingVolumetricCostEstimator(50, 0.5, 1e+4);
  AbstractCostEstimator *costEstimator = FSprout;
  GeneratorData *generatorData = new GeneratorData(16000, // Levels for tree scaling for each new segment test.
						   nFail, // Number of trials before diminish dlim.
						   lLimFactor, // Factor by which the Dlim constraint diminish after N failed trials.
						   perfusionAreaFactor, // Factor that scales the perfusion area by which Dlim is computed.
						   closeNeighborhoodFactor, // Factor that increase the neighborhood to search nearest neighbors.
						   0.1, // Factor to scale the dLim to the middle point of the new vessel to avoid close neighbors.
						   DeltaNu, // Number of bifurcation sites tested in the optimization process is given by nBifurcationTest * ( nBifurcationTest - 1 ). (default 8)
						   0,	   // Functionality of the vessel generated, important for Object trees.
						   false, // Indicates if dLimCorrectionFactor must be resetted to 1 when the stage begins.
						   costEstimator); // Cost estimator for the given stage
  DomainNVR *domain = new DomainNVR(hullVTKFilename, NVRVTKFilenames,
				    nDraw, seed, generatorData);
  domain->setIsConvexDomain(true);
  domain->setIsBifPlaneContrained(false);
  domain->setMinBifurcationAngle(0.0);
  stagedDomain->addStage(nTerms[0]+1, domain);

  cout << "Staged domain initiated." << endl;

  // Checking that the root tree's .cco file exists
  ifstream is_root_tree_correct {rootTreeFilename};
  if (!is_root_tree_correct){
    cerr << "Error: file could not be opened. The root tree's .cco file could not be found." << endl;
    exit(1);
  }
  is_root_tree_correct.close();
  
  SingleVesselCCOOTree *roottree = new SingleVesselCCOOTree(rootTreeFilename, generatorData,
							    gammas[0], deltas[0], etas[0]
							    );

  // SingleVesselCCOOTree *tree = new SingleVesselCCOOTree(x0, r0, q0,
  // 							gammas[0], deltas[0], etas[0],
  // 							p0, 1e-5,
  // 							{new GeneratorData()}
  // 							);
  
  
  VTKObjectTreeNodalWriter *treeWriter = new VTKObjectTreeNodalWriter();
  roottree->print();
  roottree->save("RootTree.cco");
  treeWriter->write("RootTree.vtp", roottree);
  
  roottree->setIsInCm(true);
  roottree->setCurrentStage(0);

  int currentStage{roottree->getCurrentStage()};
  cout << "Root tree successfully loaded... ";
  cout << "Current stage is: " << currentStage << endl; 
  
  long long int nTermTotal = roottree->getNTerms();
  cout << "Found " << nTermTotal << " terminal vessels in the root tree." << endl;
  for (int nTerm: nTerms)
    nTermTotal += nTerm;
  
  StagedFRROTreeGenerator *tree_generator = new StagedFRROTreeGenerator(stagedDomain, roottree,
									nTermTotal,
									gammas,
									deltas,
									etas);
  // tree_generator->setDLim(tree_generator->getDLim()/2.);
  
  cout << "Staged tree generator initialised." << endl;
  
  cout << "Starting tree generation." << endl;
  SingleVesselCCOOTree *tree = static_cast<SingleVesselCCOOTree *>(tree_generator->resume(200, "./"));
  cout << "Finished generating the tree." << endl;

  cout << "Saving the results..." << endl;

  tree->save(outputFilename + ".cco");  
  treeWriter->write(outputFilename + ".vtp", tree);
  cout << "Output written in " << outputFilename << ".cco and " << outputFilename << ".vtp." << endl;
}


int main(int argc, char *argv[])
{
#pragma omp parallel
  {
    cout << omp_get_num_threads() << " threads running." << endl;
  }


  vector<long long int> nTerms = {1000};
  vector<AbstractConstraintFunction<double, int>*> gammas = {new ConstantConstraintFunction<double, int>(3.0)},
    deltas = {new ConstantConstraintFunction<double, int>(0.8)},
    etas = {new ConstantConstraintFunction<double, int>(3.6)};
  double thetaMin {0}, // In radian?
    perfusionAreaFactor {.9},
    closeNeighborhoodFactor {50.0},
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
	      hullVTKFilename,
	      NVRVTKFilenames,
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
