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
#include<constrains/ConstantPiecewiseConstraintFunction.h>
#include<structures/tree/SingleVesselCCOOTree.h>
#include<structures/tree/AbstractObjectCCOTree.h>
#include<structures/vascularElements/AbstractVascularElement.h>
#include<io/VTKObjectTreeNodalWriter.h>
#include<structures/domain/DomainNVR.h>
#include<structures/domain/SimpleDomain.h>
#include<structures/domain/NormalDistributionGenerator.h>
#include<structures/vascularElements/SingleVessel.h>
#include<structures/tree/VolumetricCostEstimator.h>
#include<structures/tree/AdimSproutingVolumetricCostEstimator.h>
#include<structures/tree/SproutingVolumetricCostEstimator.h>
#include<creators/ParallelepipedCreator.h>
#include<creators/CylinderCreator.h>
#include<structures/domain/AnnulusDistributionGenerator.h>
#include<structures/CCOCommonStructures.h>

using namespace std;
// string rootTreeFilename = "CRA.cco";
// double q0 {15.0 * 0.001/60.0};	// Converts muL/min to cm^3/s
// double r0 {0.007};		// mm
// point x0 {0.26, 0.0, 0.0};	// mm
// double p0 {50 * 133.3224};  // mmHg*133.3224 = Pa, should be consistent with pressure unit 

// string outputFilename = "../Results/2DRetina",
//   hullVTKFilename = "../vtkFiles/Results/2D/Hull.vtp";
// vector<string> NVRVTKFilenames({"../vtkFiles/Results/2D/FAZ.vtp"});


void Vascularise(string outputFilename,
		 string rootTreeFilename,
		 vector<GeneratorData*> generatorData,
		 vector<long long int> nTerms,
		 vector<AbstractConstraintFunction<double, int>*> gammas,
		 vector<AbstractConstraintFunction<double, int>*> deltas,
		 vector<AbstractConstraintFunction<double, int>*> etas,
		 int nDraw, int seed,
		 double thetaMin)
{
  StagedDomain *stagedDomain = new StagedDomain();

  // Checking that the root tree's .cco file exists
  ifstream is_root_tree_correct {rootTreeFilename};
  if (!is_root_tree_correct){
    cerr << "Error: file could not be opened. The root tree's .cco file could not be found." << endl;
    exit(1);
  }
  is_root_tree_correct.close();
  
  SingleVesselCCOOTree *rootTree = new SingleVesselCCOOTree(rootTreeFilename, generatorData[0],
							    gammas[0], deltas[0], etas[0]
							    );


  // SCP - Whole Retina
  
  vector<vector<double>> vertices = rootTree->getVertices();
  double lb[3] {-2.0,-2.0,-0.01}, ub[3] {2.0,2.0,0.0};
  
  point opticDisc = rootTree->getRoot()->getVessels()[0]->xDist;
  cout << opticDisc << endl;
  SCPDistributionGenerator *dist = new SCPDistributionGenerator(0.05, opticDisc.p, -0.5,4, -0.008, 0.00);

  ParallelepipedCreator *para = new ParallelepipedCreator(lb, ub);
  para->create("tmp.vtk");
  SimpleDomain *domainSCP = new SimpleDomain("tmp.vtk", nDraw, seed,
					  generatorData[0], dist);
  domainSCP->setIsConvexDomain(true);
  domainSCP->setIsBifPlaneContrained(false);
  domainSCP->setMinBifurcationAngle(thetaMin);
  domainSCP->setOpticDisc(opticDisc);
  stagedDomain->addStage(nTerms[0]+1, domainSCP);
  stagedDomain->setOpticDisc(opticDisc);  


  // SCP - MACULA

  delete para;
  // Change the FOV
  for (int i = 0; i<2; i++){
    lb[i] = -0.4;
    ub[i] = 0.4;
  }  
  para = new ParallelepipedCreator(lb, ub);
  para->create("tmp.vtk");

  SCPDistributionGenerator *distMacula = new SCPDistributionGenerator(0.05, opticDisc.p, -0.01,4, -0.008, 0.00);
  SimpleDomain *domainSCPMacula = new SimpleDomain("tmp.vtk", nDraw, seed,
						   generatorData[1], distMacula);
  domainSCPMacula->setIsConvexDomain(true);
  domainSCPMacula->setIsBifPlaneContrained(false);
  domainSCPMacula->setMinBifurcationAngle(thetaMin);
  domainSCPMacula->setOpticDisc(opticDisc);
  stagedDomain->addStage(nTerms[1]+1, domainSCPMacula);  
  
  cout << "Staged domain initiated." << endl;


  VTKObjectTreeNodalWriter *treeWriter = new VTKObjectTreeNodalWriter();
  rootTree->print();
  rootTree->save("RootTree.cco");
  treeWriter->write("RootTree.vtp", rootTree);
  
  rootTree->setIsInCm(true);
  rootTree->setCurrentStage(1);

  int currentStage{rootTree->getCurrentStage()};
  cout << "Root tree successfully loaded... ";
  cout << "Current stage is: " << currentStage << endl; 
  
  long long int nTermTotal = rootTree->getNTerms();
  cout << "Found " << nTermTotal << " terminal vessels in the root tree." << endl;
  for (int nTerm: nTerms)
    nTermTotal += nTerm;
  
  StagedFRROTreeGenerator *tree_generator = new StagedFRROTreeGenerator(stagedDomain, rootTree,
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
  omp_set_dynamic(0);
  omp_set_num_threads(24);
  cout << omp_get_num_threads() << " threads running." << endl;
  cout << omp_get_max_threads() << " threads available." << endl;

  string outputFilename {"../Results/3DSCP_AllAtOnce"},
    inputCCO {"./base_trees/SegmentedArcades_12.cco"};


  // VolumetricCostEstimator *FSprout = new VolumetricCostEstimator();
  // SproutingVolumetricCostEstimator *FSprout = new SproutingVolumetricCostEstimator(50, 0.5, 1e+4);
  AdimSproutingVolumetricCostEstimator *FSprout1 = new AdimSproutingVolumetricCostEstimator(50, 0.5, 1e+4, 40, 70*1e-6);
  AdimSproutingVolumetricCostEstimator *FSprout2 = new AdimSproutingVolumetricCostEstimator(50, 0.5, 1e+4, 40, 70*1e-6);

  AbstractCostEstimator *costEstimator1 = FSprout1, *costEstimator2 = FSprout2;
  
  GeneratorData *generatorData1 = new GeneratorData(16000, // Levels for tree scaling for each new segment test.
						    20, // Number of trials before diminish dlim.
						    0.5, // Factor by which the Dlim constraint diminish after N failed trials.
						    3.0, // Factor that scales the perfusion area by which Dlim is computed.
						    1.0, // Factor that increase the neighborhood to search nearest neighbors.
						    0.1, // Factor to scale the dLim to the middle point of the new vessel to avoid close neighbors.
						    8, // Number of bifurcation sites tested in the optimization process is given by nBifurcationTest * ( nBifurcationTest - 1 ). (default 8)
						    0,	   // Functionality of the vessel generated, important for Object trees.
						    false, // Indicates if dLimCorrectionFactor must be resetted to 1 when the stage begins.
						    costEstimator1); // Cost estimator for the given stage

  GeneratorData *generatorData2 = new GeneratorData(16000, // Levels for tree scaling for each new segment test.
						    20, // Number of trials before diminish dlim.
						    0.5, // Factor by which the Dlim constraint diminish after N failed trials.
						    3.0, // Factor that scales the perfusion area by which Dlim is computed.
						    1.0, // Factor that increase the neighborhood to search nearest neighbors.
						    0.1, // Factor to scale the dLim to the middle point of the new vessel to avoid close neighbors.
						    8, // Number of bifurcation sites tested in the optimization process is given by nBifurcationTest * ( nBifurcationTest - 1 ). (default 8)
						    0,	   // Functionality of the vessel generated, important for Object trees.
						    false, // Indicates if dLimCorrectionFactor must be resetted to 1 when the stage begins.
						    costEstimator2); // Cost estimator for the given stage

  vector<GeneratorData*> generatorData = {generatorData1, generatorData2};
  
  // Consecutive attempts to generate a point - nFail
  int nFail = 20;
  // Discretisation of the testing triangle for the bifurcation - Delta nu - Figure 1
  int DeltaNu = 7;
  // Buffer size for random point generation
  int nDraw {1000};
  // Random seed
  long long int seed {time(nullptr)};

  AbstractConstraintFunction<double, int> *gamma1 {new ConstantConstraintFunction<double, int>(2.85)}; // Murray's law exponent
  AbstractConstraintFunction<double, int> *gamma2 {new ConstantConstraintFunction<double, int>(3.0)}; // Murray's law exponent

  AbstractConstraintFunction<double, int> *delta1 {new ConstantPiecewiseConstraintFunction<double, int>({0.1, 0.8, 0.4}, {0, 10, 100})};
  AbstractConstraintFunction<double, int> *delta2 {new ConstantPiecewiseConstraintFunction<double, int>({0.1, 0.8, 0.4}, {0, 10, 100})};
  
  AbstractConstraintFunction<double, int> *eta1 {new ConstantConstraintFunction<double, int>(0.36)}; // viscosity
  AbstractConstraintFunction<double, int> *eta2 {new ConstantConstraintFunction<double, int>(0.5)}; // viscosity

  vector<AbstractConstraintFunction<double, int>*> gammas = {gamma1,gamma2},
    deltas = {delta1, delta2},
    etas = {eta1, eta2};

  double thetaMin = 0.16 * M_PI;

  
  Vascularise(outputFilename,
	      inputCCO,
	      generatorData, // GeneratorData
	      {250,250}, // nTerms
	      gammas,
	      deltas,
	      etas,
	      10000, seed,
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
