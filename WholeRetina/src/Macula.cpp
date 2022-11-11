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
#include<structures/domain/PartiallyVascularizedDomain.h>
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
		 string Omega1,
		 string Omega2,
		 string FAZ, 
		 std::vector<long long int> nTerms,
		 std::vector<AbstractConstraintFunction<double, int>*> gammas,
		 std::vector<AbstractConstraintFunction<double, int>*> deltas,
		 std::vector<AbstractConstraintFunction<double, int>*> etas,
		 int nDraw, std::vector<long long int> seeds, int nFail, double lLimFactor,
		 double perfusionAreaFactor, double closeNeighborhoodFactor, int DeltaNu,
		 std::vector<double> thetaMin)
{  

  std::cout << "The .vtp files: " << hullVTKFilename << " " << Omega1 << " " << Omega2 << " " << FAZ << std::endl;
  
  // VolumetricCostEstimator *FSprout = new VolumetricCostEstimator();
  SproutingVolumetricCostEstimator *FSprout = new SproutingVolumetricCostEstimator(10, 0.5, 1e+4);
  AbstractCostEstimator *costEstimator = FSprout;

  // This generator is only for importing the root tree and is not used in the generation.
  GeneratorData *generatorData = new GeneratorData(16000, nFail,lLimFactor, perfusionAreaFactor, closeNeighborhoodFactor, 0.1, DeltaNu, 0, false, costEstimator);
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
  VTKObjectTreeNodalWriter *treeWriter = new VTKObjectTreeNodalWriter();
  treeWriter->write("RootTreeMacula.vtp", rootTree);
  rootTree->save("RootTreeMacula.cco");
  rootTree->setIsInCm(true);
  int nTermRoot = rootTree->getNTerms();
  int stageRoot = rootTree->getCurrentStage();
  std::cout << "The root tree has " << nTermRoot << " terminals and current stage is " << stageRoot << std::endl;
  


  
  // // // Stage 1, within an annulus in the outer macula
  GeneratorData *generatorData1 = new GeneratorData(16000, // Levels for tree scaling for each new segment test.
						    20, // Number of trials before diminish dlim.
						    0.9, // Factor by which the Dlim constraint diminish after N failed trials.
						    perfusionAreaFactor, // Factor that scales the perfusion area by which Dlim is computed.
						    closeNeighborhoodFactor, // Factor that increase the neighborhood to search nearest neighbors.
						    0.1, // Factor to scale the dLim to the middle point of the new vessel to avoid close neighbors.
						    DeltaNu, // Number of bifurcation sites tested in the optimization process is given by nBifurcationTest * ( nBifurcationTest - 1 ).
						    0,	   // Functionality of the vessel generated, important for Object trees.
						    false, // Indicates if dLimCorrectionFactor must be resetted to 1 when the stage begins.
						    costEstimator); // Cost estimator for the given stage
  // The domain for stage 1
  std::vector<string> VRVTK1 = {Omega1},
    NVRVTK1 = {Omega2, FAZ};

  
  PartiallyVascularizedDomain *domain1 = new PartiallyVascularizedDomain(hullVTKFilename, VRVTK1, NVRVTK1,
									 nDraw, seeds[0], generatorData1);
  domain1->setIsConvexDomain(true);
  domain1->setIsBifPlaneContrained(false);
  domain1->setMinBifurcationAngle(thetaMin[0]);


  int nTermTotal = nTermRoot + nTerms[0];
  StagedDomain *stagedDomain = new StagedDomain();
  stagedDomain->setInitialStage(stageRoot+1);
  
  stagedDomain->addStage(nTerms[0], domain1);

  StagedFRROTreeGenerator *treeGenerator = new StagedFRROTreeGenerator(stagedDomain,
								       rootTree,
								       nTermTotal,
								       {gammas[0]},
								       {deltas[0]},
								       {etas[0]});
  // // Generate the new macular vessels

  rootTree = static_cast<SingleVesselCCOOTree *>(treeGenerator->resume(10, "./LogCCO/"));
  std::cout << "Finished generating stage 1." << std::endl;
  
  // // Save
  std::cout << "Saving the results." << std::endl;
  rootTree->save(outputFilename + "_root_with_macula_stage_1.cco");
  treeWriter->write(outputFilename + "_root_with_macula_stage_1.vtp", rootTree);

  delete treeGenerator;
  delete stagedDomain;
  
  
  // // // Stage 2, within the macula but outside the annulus of stage 1
  GeneratorData *generatorData2 = new GeneratorData(16000, // Levels for tree scaling for each new segment test.
						    nFail, // Number of trials before diminish dlim.
						    lLimFactor, // Factor by which the Dlim constraint diminish after N failed trials.
						    perfusionAreaFactor, // Factor that scales the perfusion area by which Dlim is computed.
						    closeNeighborhoodFactor, // Factor that increase the neighborhood to search nearest neighbors.
						    0.1, // Factor to scale the dLim to the middle point of the new vessel to avoid close neighbors.
						    DeltaNu, // Number of bifurcation sites tested in the optimization process is given by nBifurcationTest * ( nBifurcationTest - 1 ).
						    0,	   // Functionality of the vessel generated, important for Object trees.
						    false, // Indicates if dLimCorrectionFactor must be resetted to 1 when the stage begins.
						    costEstimator); // Cost estimator for the given stage
  // The domain for stage 2
  std::vector<string> VRVTK2 = {Omega2},
    NVRVTK2 = {Omega1, FAZ};

  PartiallyVascularizedDomain *domain2 = new PartiallyVascularizedDomain(hullVTKFilename, VRVTK2, NVRVTK2,
									 nDraw, seeds[1], generatorData2);
  domain2->setIsConvexDomain(true);
  domain2->setIsBifPlaneContrained(false);
  domain2->setMinBifurcationAngle(thetaMin[1]);
 
  nTermTotal = rootTree->getNTerms() + nTerms[1];
  stageRoot  = rootTree->getCurrentStage();
  stagedDomain = new StagedDomain();
  stagedDomain->setInitialStage(stageRoot+1); 
  stagedDomain->addStage(nTerms[1], domain2);

  treeGenerator = new StagedFRROTreeGenerator(stagedDomain,
					      rootTree,
					      nTermTotal,
					      {gammas[1]},
					      {deltas[1]},
					      {etas[1]});
  // // Generate the new macular vessels
  rootTree = static_cast<SingleVesselCCOOTree *>(treeGenerator->resume(10, "./LogCCO/"));
  std::cout << "Finished generating stage 2." << std::endl;
  
  // // Save
  std::cout << "Saving the results." << std::endl;
  rootTree->save(outputFilename + "_root_with_macula_stage_2.cco");
  treeWriter->write(outputFilename + "_root_with_macula_stage_2.vtp", rootTree);

  delete treeGenerator;
  delete stagedDomain;

  

  // // // Stage 3, within the macula but outside the annulus of stage 1
  GeneratorData *generatorData3 = new GeneratorData(16000, // Levels for tree scaling for each new segment test.
						    nFail/5, // Number of trials before diminish dlim.
						    lLimFactor, // Factor by which the Dlim constraint diminish after N failed trials.
						    perfusionAreaFactor, // Factor that scales the perfusion area by which Dlim is computed.
						    closeNeighborhoodFactor, // Factor that increase the neighborhood to search nearest neighbors.
						    0.1, // Factor to scale the dLim to the middle point of the new vessel to avoid close neighbors.
						    DeltaNu, // Number of bifurcation sites tested in the optimization process is given by nBifurcationTest * ( nBifurcationTest - 1 ).
						    0,	   // Functionality of the vessel generated, important for Object trees.
						    false, // Indicates if dLimCorrectionFactor must be resetted to 1 when the stage begins.
						    costEstimator); // Cost estimator for the given stage
  // The domain for stage 3
  std::vector<string> VRVTK3 = {Omega1, Omega2},
    NVRVTK3 = {FAZ};
  
  PartiallyVascularizedDomain *domain3 = new PartiallyVascularizedDomain(hullVTKFilename, VRVTK3, NVRVTK3,
									 nDraw, seeds[2], generatorData3);
  domain3->setIsConvexDomain(true);
  domain3->setIsBifPlaneContrained(false);
  domain3->setMinBifurcationAngle(thetaMin[2]);  

  nTermTotal = rootTree->getNTerms() + nTerms[2];
  stageRoot  = rootTree->getCurrentStage();
  stagedDomain = new StagedDomain();
  stagedDomain->setInitialStage(stageRoot+1); 
  stagedDomain->addStage(nTerms[2], domain3);

  treeGenerator = new StagedFRROTreeGenerator(stagedDomain,
					      rootTree,
					      nTermTotal,
					      {gammas[2]},
					      {deltas[2]},
					      {etas[2]});
  // // Generate the new macular vessels

  rootTree = static_cast<SingleVesselCCOOTree *>(treeGenerator->resume(10, "./LogCCO/"));
  std::cout << "Finished generating stage 3." << std::endl;
  
  // // Save
  std::cout << "Saving the results." << std::endl;
  rootTree->save(outputFilename + "_root_with_macula_stage_3.cco");
  treeWriter->write(outputFilename + "_root_with_macula_stage_3.vtp", rootTree);

  delete treeGenerator;
  delete stagedDomain;

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
  std::vector<string> NVRVTKFilenames;
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
  AbstractConstraintFunction<double, int> *delta0 {new ConstantConstraintFunction<double, int>(stod(line))}; // Symmetry ratio
  // AbstractConstraintFunction<double, int> *delta0 {new ConstantPiecewiseConstraintFunction<double, int>({0.9, 0.}, {0, 5})};
  AbstractConstraintFunction<double, int> *delta1 {new ConstantPiecewiseConstraintFunction<double, int>({0.9, 0.5}, {0, 5})};
  AbstractConstraintFunction<double, int> *delta2 {new ConstantPiecewiseConstraintFunction<double, int>({0.9, 0.5}, {0, 5})};
  
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  AbstractConstraintFunction<double, int> *eta {new ConstantConstraintFunction<double, int>(stod(line))}; // viscosity

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double viscTolerance {stod(line)};

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double fractionOfPi {stod(line)};
  std::vector<double> thetaMin = {fractionOfPi * M_PI, fractionOfPi * M_PI, 2*fractionOfPi * M_PI};
  
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
  
  std::vector<AbstractConstraintFunction<double, int>*> gammas = {gamma, gamma, gamma},
    deltas = {delta0, delta1, delta2},
    etas = {eta, eta, eta};
  
  // Consecutive attempts to generate a point - nFail
  int nFail = 200;
  // Discretisation of the testing triangle for the bifurcation - Delta nu - Figure 1
  int DeltaNu = 7;
  // Buffer size for random point generation
  int nDraw {1000};
  // Random seed
  std::vector<long long int> seeds {time(nullptr), time(nullptr), time(nullptr)};
  std::cout << "The random seeds are " << seeds[0] << ", " << seeds[1] << " and " << seeds[2] << std::endl;


  string Omega1 = "../vtkFiles/Results/2D/Macula_omega1.vtp",
    Omega2 = "../vtkFiles/Results/2D/Macula_omega2.vtp",
    FAZ = "../vtkFiles/Results/2D/FAZ.vtp";
  
  Vascularise(outputFilename,
	      inputCCO,
	      hullVTKFilename,
	      Omega1,
	      Omega2,
	      FAZ,
	      {20, 30, 200},
	      gammas,
	      deltas,
	      etas,
	      nDraw, seeds, nFail, lLimFactor,
	      perfusionAreaFactor, closeNeighborhoodFactor, DeltaNu,
	      thetaMin);

  // for (auto p : gammas)
  //   delete p;
  // gammas.clear();
  // for (auto p : etas)
  //   delete p;
  // etas.clear();
  // for (auto p : deltas)
  //   delete p;
  // deltas.clear();
  return 0;
} 
