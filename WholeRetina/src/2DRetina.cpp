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
#include<structures/domain/SimpleDomain2D.h>
#include<structures/domain/NormalDistributionGenerator.h>
#include<structures/vascularElements/SingleVessel.h>
#include<structures/tree/VolumetricCostEstimator.h>
#include<structures/tree/AdimSproutingVolumetricCostEstimator.h>
#include<structures/tree/SproutingVolumetricCostEstimator.h>
#include<creators/ParallelepipedCreator.h>
#include<creators/CylinderCreator.h>

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




  VolumetricCostEstimator *FSprout = new VolumetricCostEstimator();
  // SproutingVolumetricCostEstimator *FSprout = new SproutingVolumetricCostEstimator(50, 0.5, 1e+4);
  // AdimSproutingVolumetricCostEstimator *FSprout = new AdimSproutingVolumetricCostEstimator(50, 0.5, 1e+4, 40, 70*1e-6);
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
  domain->setMinBifurcationAngle(thetaMin);
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
  cout << omp_get_num_threads() << " threads running." << endl;

  // Read the configuration file passed as command line argument
  string ConfigurationFileName = argv[1];
  ifstream config;
  config.open(ConfigurationFileName);
  if (!config){
    cerr << "Error: parameter file " << ConfigurationFileName
	 << " could not be opened." << endl;
    exit(1);
  }
  else
    cout << "Reading parameters from " << ConfigurationFileName << endl;

  // Reading input and output files and parameters
  
  // Output file
  string line;
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  string outputFilename {line};
  
  // Hull
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  string hullVTKFilename {line};
  cout << "Hull file: " << hullVTKFilename << endl;				       

  // Domain 1
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  int nNVR {stoi(line)};  
  vector<string> NVRVTKFilenames;
  cout << "The " << nNVR << " non vascular regions for domain 1 are: ";
  for (int i = 0; i < nNVR; i++){
    getline(config, line);
    NVRVTKFilenames.push_back(line);
    cout << NVRVTKFilenames[i] << " ";
  }
  cout << endl;

  // Root tree
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  string inputCCO {line};
  cout << "Root tree file: " << inputCCO << endl;
  
  // Other simulation parameters
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  long long int nTerms {stoll(line)};
  
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double lLimFactor {stod(line)};

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  AbstractConstraintFunction<double, int> *gamma {new ConstantConstraintFunction<double, int>(stod(line))}; // Murray's law exponent

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  // AbstractConstraintFunction<double, int> *delta {new ConstantConstraintFunction<double, int>(stod(line))}; // Symmetry ratio
  AbstractConstraintFunction<double, int> *delta {new ConstantPiecewiseConstraintFunction<double, int>({0.9, 0.8, 0.4}, {0, 2, 5})};
  
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  AbstractConstraintFunction<double, int> *eta {new ConstantConstraintFunction<double, int>(stod(line))}; // viscosity

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double viscTolerance {stod(line)};

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double fractionOfPi {stod(line)};
  double thetaMin = fractionOfPi * M_PI;
  
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double perfusionAreaFactor {stod(line)};

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double closeNeighborhoodFactor {stod(line)};

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  double q0, p0, rEyeball;
  config >> q0 >> p0;
  config >> rEyeball;
  cout << "Inlet flow and pressure: " << q0 << " " << p0 << " in a retina of radius " << rEyeball << "cm" << endl;
  
  
  // Copying the config file in the output folder
  string root_results {outputFilename.substr(0, outputFilename.find_last_of("/"))};
  string filename {ConfigurationFileName.substr(ConfigurationFileName.find_last_of("/") + 1)};
  string command {"mkdir -p " + root_results};
  const int dir_err = system(command.c_str());
  if (-1 == dir_err)
    {
      printf("Error creating directory!\n");
      exit(1);
    }
  command = {"cp " + ConfigurationFileName + " " + root_results + "/" + filename};
  const int copy_err = system(command.c_str());
  if (-1 == copy_err)
    {
      printf("Error copying the config file!\n");
      exit(1);
    }
  
  cout << "Simulation parameters read successfully." << endl;
  
  vector<AbstractConstraintFunction<double, int>*> gammas = {gamma},
    deltas = {delta},
    etas = {eta};
  
  // Consecutive attempts to generate a point - nFail
  int nFail = 200;
  // Discretisation of the testing triangle for the bifurcation - Delta nu - Figure 1
  int DeltaNu = 7;
  // Buffer size for random point generation
  int nDraw {10000};
  // Random seed
  long long int seed {time(nullptr)};

  double lb[3] = {-1., -1., 0.0}, ub[3] = {1.0,1.0,0.005};
  ParallelepipedCreator *Hull = new ParallelepipedCreator(lb, ub);
  Hull->create(hullVTKFilename);
  vector<double> center = {0.188*1.5, 0.0, 0.0};
  double radiusMacula   = 0.3 / 2.0;
  CylinderCreator *Macula = new CylinderCreator(center, radiusMacula, 0.005, 10);
  Macula->create(NVRVTKFilenames[0]);
  
  Vascularise(outputFilename,
	      inputCCO,
	      hullVTKFilename,
	      NVRVTKFilenames,
	      {nTerms},
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
