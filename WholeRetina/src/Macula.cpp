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


  // //    Import the root tree
  // Checking that the root tree's .cco file exists
  ifstream is_root_tree_correct {rootTreeFilename};
  if (!is_root_tree_correct){
    cerr << "Error: file could not be opened. The root tree's .cco file could not be found." << endl;
    exit(1);
  }
  is_root_tree_correct.close();
  std::cout << "Root tree file found, reading the root tree." << std::endl;
  SingleVesselCCOOTree *rootTree = new SingleVesselCCOOTree(rootTreeFilename, generatorData,
							    gammas[0], deltas[0], etas[0]
							    );
  rootTree->setIsInCm(true);
  int nTermRoot = rootTree->getNTerms();
  int stageRoot = rootTree->getCurrentStage();
  std::cout << "The root tree has " << nTermRoot << " terminals and current stage is " << stageRoot << std::endl;
  
  // // The domain and generator
  StagedDomain *stagedDomain = new StagedDomain();
  DomainNVR *domain = new DomainNVR(hullVTKFilename, NVRVTKFilenames,
				    nDraw, seed, generatorData);
  domain->setIsConvexDomain(true);
  domain->setIsBifPlaneContrained(false);
  domain->setMinBifurcationAngle(thetaMin);
  stagedDomain->setInitialStage(stageRoot+1);
  stagedDomain->addStage(nTerms[0] + 1, domain);

  int nTermTotal = nTermRoot + nTerms[0];
  
  StagedFRROTreeGenerator *treeGenerator = new StagedFRROTreeGenerator(stagedDomain,
								       rootTree,
								       nTermTotal,
								       gammas,
								       deltas,
								       etas);
  std::cout << "Initial DLim = " << treeGenerator->getDLim() << std::endl;
  treeGenerator->setDLim(treeGenerator->getDLim()/4.0);
  std::cout << "Ready to generate the macular network." << std::endl;

  // // Generate the new macular vessels

  SingleVesselCCOOTree *tree = static_cast<SingleVesselCCOOTree *>(treeGenerator->resume(200, "./"));
  std::cout << "Finished generating the tree." << std::endl;

  // // Save
  
  VTKObjectTreeNodalWriter *treeWriter = new VTKObjectTreeNodalWriter();
  std::cout << "Printing the root tree." << std::endl;
  rootTree->print();

  std::cout << "Saving the results." << std::endl;
  tree->save(outputFilename + "_root_with_macula.cco");
  treeWriter->write(outputFilename + "_root_with_macula.vtp", tree);

  std::cout << "Output writen in " << outputFilename + "_root_with_macula.vtp/cco" << std::endl;
  
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
  AbstractConstraintFunction<double, int> *delta {new ConstantPiecewiseConstraintFunction<double, int>({0.4, 0.8}, {0, 5})};
  
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
