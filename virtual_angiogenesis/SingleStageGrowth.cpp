// My libs
#include<IntercapillaryDistance.hpp>
#include<VascularIndexes.hpp>

// STD Libs
#include<cmath>
#include<string>
#include<ctime>
#include<cstdio>
#include<cstdlib>
#include<vector>
#include<fstream>
#include<omp.h>

// vtk libs
#include<vtkSmartPointer.h>
#include<vtkPolyData.h>

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
#include<structures/tree/VolumetricCostEstimator.h>

// For tree stats
#include<stats/ObjectTreeStatsManager.h>
#include<stats/VesselObjectHandler.h>

using namespace std;

void Vascularise(string output_filename, string rootTreeFilename, string Hull,
		 vector<string> NVR,
		 long long int nTerm,
		 AbstractConstraintFunction<double, int>* gam,
		 AbstractConstraintFunction<double, int>* delta,
		 AbstractConstraintFunction<double, int>* eta,
		 int nDraw, int seed, int nFails, double lLimFactor,
		 double perfusion_area_factor, double closeNeighborhoodFactor, int nu,
		 double theta_min,
		 double targetICD)
{
  
  // SproutingVolumetricCostEstimator *FSprout = new SproutingVolumetricCostEstimator(50, 0.5, 1e+4);
  VolumetricCostEstimator *FSprout = new VolumetricCostEstimator();
  AbstractCostEstimator *costEstimator = FSprout;
  GeneratorData *genData = new GeneratorData(160, nFails, lLimFactor, perfusion_area_factor,
					      closeNeighborhoodFactor, 0.25, nu, 0, false,
					      costEstimator);
  
  // Domain definition for stage 
  DomainNVR *domain = new DomainNVR(Hull, NVR, nDraw, seed, genData);
  domain->setIsConvexDomain(false);
  domain->setMinBifurcationAngle(theta_min);
  cout << "Domain generated." << endl;

  StagedDomain *stagedDomain = new StagedDomain();
  stagedDomain->addStage(nTerm, domain);
  cout << "Staged Domain initialized.";
  
  // Checking that the root tree's .cco file exists
  ifstream isRootTreeCorrect {rootTreeFilename};
  if (!isRootTreeCorrect){
    cerr << "Error: file could not be opened. The root tree's .cco file could not be found." << endl;
    exit(1);
  }
  isRootTreeCorrect.close();
  
  SingleVesselCCOOTree *tree = new SingleVesselCCOOTree(rootTreeFilename, genData, gam, delta, eta);
  tree->setIsInCm(true);
  tree->setCurrentStage(2);
  int currentStage{tree->getCurrentStage()};
  cout << "Root tree successfully loaded... ";
  cout << "Current stage is: " << currentStage << endl; 
  cout << "Saving root tree.\n" << endl;
  // Save the root tree
  VTKObjectTreeNodalWriter *tree_writer = new VTKObjectTreeNodalWriter();
  tree_writer->write(output_filename + "_root.vtp", tree);
  tree->save(output_filename + ".cco.root");
  stagedDomain->setInitialStage(0);

  long long int nTermTotal = nTerm + tree->getNTerms();
  
  StagedFRROTreeGenerator *treeGenerator = new StagedFRROTreeGenerator(stagedDomain, tree,
								       nTermTotal,
								       {gam},
								       {delta},
								       {eta});
  treeGenerator->setDLim(treeGenerator->getDLim()/2.);  
  cout << "Staged tree generator initialized." << endl;
  
  cout << "Starting tree generation." << endl;
  tree = {(SingleVesselCCOOTree *) treeGenerator->resume(200, "./")};

  
  vtkSmartPointer<vtkPolyData> treePolyData = tree->getVtkTree();
  // bool isCriterionReached = ( targetICD>IntercapillaryDistance(treePolyData, 0.3, 1024) );
  // cout << "Intercapillary Distance  |" << IntercapillaryDistance(tree->getVtkTree(), 0.3, 512) << endl;
  // cout << "Vessel Area Density      |" << VesselAreaDensity(tree->getVessels(), 0.09-pow(0.04, 2)*M_PI) << endl;
  // cout << "Vessel Skeleton Density  |" << VesselSkeletonDensity(tree->getVessels(), 0.09-pow(0.04, 2)*M_PI) << endl;
  // cout << "Vessel Perimeter Index   |" << VesselPerimeterIndex(tree->getVessels(), 0.09-pow(0.04, 2)*M_PI) << endl;
  // cout << "Vessel Complexity Index  |" << VesselComplexityIndex(tree->getVessels(), 0.09-pow(0.04, 2)*M_PI) << endl;
  // cout << "Vessel Diameter Index    |" << VesselDiameterIndex(tree->getVessels(), 0.09-pow(0.04, 2)*M_PI) << endl;

  vector<vector<double>> metricsObserver;

  {
    vector<double> metrics;
    metrics.push_back(tree->getNTerms()); 
    metrics.push_back(treeGenerator->getDLim());
    metrics.push_back(tree->computeTreeCost(tree->getRoot()));
    metrics.push_back(( tree->getConnectivity() ).size()); 
    metrics.push_back(IntercapillaryDistance(tree->getVtkTree(), 0.3, 512));
    metrics.push_back(VesselAreaDensity(tree->getVessels(), 0.09-pow(0.04, 2)*M_PI));
    metrics.push_back(VesselSkeletonDensity(tree->getVessels(), 0.09-pow(0.04, 2)*M_PI));
    metrics.push_back(VesselPerimeterIndex(tree->getVessels(), 0.09-pow(0.04, 2)*M_PI));
    metrics.push_back(VesselComplexityIndex(tree->getVessels(), 0.09-pow(0.04, 2)*M_PI));
    metrics.push_back(VesselDiameterIndex(tree->getVessels(), 0.09-pow(0.04, 2)*M_PI));
    
    metricsObserver.push_back(metrics);
    // metrics.empty();
  }
  
  for (int i = 0; i < 10; i++)
    {
      double dLim = treeGenerator->getDLim();
      stagedDomain = new StagedDomain();
      stagedDomain->addStage(50, domain);
      nTermTotal = nTermTotal+50;
      treeGenerator = new StagedFRROTreeGenerator(stagedDomain, tree,
						  nTermTotal,
						  {gam},
						  {delta},
						  {eta});
      tree = {(SingleVesselCCOOTree *) treeGenerator->resume(200, "./")};

      // cout << "Intercapillary Distance  |" << IntercapillaryDistance(tree->getVtkTree(), 0.3, 512) << endl;
      // cout << "Vessel Area Density      |" << VesselAreaDensity(tree->getVessels(), 0.09-pow(0.04, 2)*M_PI) << endl;
      // cout << "Vessel Skeleton Density  |" << VesselSkeletonDensity(tree->getVessels(), 0.09-pow(0.04, 2)*M_PI) << endl;
      // cout << "Vessel Perimeter Index   |" << VesselPerimeterIndex(tree->getVessels(), 0.09-pow(0.04, 2)*M_PI) << endl;
      // cout << "Vessel Complexity Index  |" << VesselComplexityIndex(tree->getVessels(), 0.09-pow(0.04, 2)*M_PI) << endl;
      // cout << "Vessel Diameter Index    |" << VesselDiameterIndex(tree->getVessels(), 0.09-pow(0.04, 2)*M_PI) << endl;
      // cout << "Area of the domain       |" << 0.09-pow(0.04, 2)*M_PI << endl;

      vector<double> metrics;
      metrics.push_back(tree->getNTerms()); 
      metrics.push_back(treeGenerator->getDLim());
      metrics.push_back(tree->computeTreeCost(tree->getRoot()));
      metrics.push_back(( tree->getConnectivity() ).size()); 
      metrics.push_back(IntercapillaryDistance(tree->getVtkTree(), 0.3, 512));
      metrics.push_back(VesselAreaDensity(tree->getVessels(), 0.09-pow(0.04, 2)*M_PI));
      metrics.push_back(VesselSkeletonDensity(tree->getVessels(), 0.09-pow(0.04, 2)*M_PI));
      metrics.push_back(VesselPerimeterIndex(tree->getVessels(), 0.09-pow(0.04, 2)*M_PI));
      metrics.push_back(VesselComplexityIndex(tree->getVessels(), 0.09-pow(0.04, 2)*M_PI));
      metrics.push_back(VesselDiameterIndex(tree->getVessels(), 0.09-pow(0.04, 2)*M_PI));
      
      metricsObserver.push_back(metrics);
      // metrics.clear();
    }

  cout << "Saving the results..." << endl;

  tree->save(output_filename + ".cco");  
  tree_writer->write(output_filename + ".vtp", tree);
  cout << "Output written in " << output_filename << ".cco and " << output_filename << ".vtp." << endl;

  ObjectTreeStatsManager *statsManager = new ObjectTreeStatsManager(tree);
  vector<double> levels, means, stds;
  VesselObjectHandler::ATTRIBUTE att = VesselObjectHandler::ATTRIBUTE::DIAMETER;
  
  // Load the vectors with levels and the diameter for each vessel
  statsManager->getMeanPerLevel(&levels, &means, &stds, att);
  string fileName = output_filename + "_diameters.dat";
  ofstream outputFile(fileName);
  for (int i = 0; i<levels.size(); i++)
    {
      outputFile << levels[i] << ' ' << means[i] << ' ' << stds[i] << endl;
    }
  outputFile.close();

  // Save the metrics at each stage
  fileName = output_filename + "_metrics.dat";
  outputFile.open(fileName);
  outputFile << "# NTerms DLim TreeCost NVessels ICD VAD VSD VPI VCI VDI" << endl;
  for (auto e: metricsObserver)
    {
      for (int j = 0; j < e.size(); j++)
	outputFile << e[j] << ' ';
      outputFile << endl;
    }
  outputFile.close();
}

int main(int argc, char *argv[])
{
#pragma omp parallel
  {
    cout << "Thread " << omp_get_thread_num() << " running." << endl;
  }
  
  // Read the configuration file passed as command line argument
  string ConfigurationFileName = argv[1];
  ifstream config;
  config.open(ConfigurationFileName);
  if (!config){
    cerr << "Error: file could not be opened." << endl;
    exit(1);
  }

  // Reading input files and parameters - output file
  string line;
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  string output_filename {line};
  
  // Hull (macula)
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  string Hull {line};
  cout << "Hull file: " << Hull << endl;				       

  // Domain 
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  int nbNVR {stoi(line)};  
  vector<string> NVR;
  cout << "The " << nbNVR << " non vascular regions for domain 1 are: ";
  for (int i = 0; i < nbNVR; i++){
    getline(config, line);
    NVR.push_back(line);
    cout << NVR[i] << " ";
  }
  cout << endl;

  // Root tree
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  string input_cco {line};
  cout << "Root tree file: " << input_cco << endl;
  
  // Other simulation parameters
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  long long int nTerm {stoll(line)};
  
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double lLimFactor {stod(line)};

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  AbstractConstraintFunction<double, int> *gam {new ConstantConstraintFunction<double, int>(stod(line))};

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  AbstractConstraintFunction<double, int> *delta {new ConstantConstraintFunction<double, int>(stod(line))};

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  AbstractConstraintFunction<double, int> *eta {new ConstantConstraintFunction<double, int>(stod(line))};

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double viscTolerance {stod(line)};

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double fraction_of_pi {stod(line)};
  double theta_min = fraction_of_pi * M_PI;
  
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double perfusion_area_factor {stod(line)};

  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double closeNeighborhoodFactor {stod(line)};

  // Intercapillary distance target
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double targetICD {stod(line)};

  
  // Copying the config file in the output folder
  string root_results {output_filename.substr(0, output_filename.find_last_of("/"))};
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

  // Consecutive attempts to generate a point - nFails
  int nFails = 200;
  // Discretisation of the testing triangle for the bifurcation - Delta nu - Figure 1
  int nu = 7;
  // Buffer size for random point generation
  int nDraw {10000};
  // Random seed
  long long int seed {time(nullptr)};
  
  Vascularise(output_filename, input_cco, Hull, NVR, nTerm,
	      gam, delta, eta, nDraw, seed, nFails, lLimFactor,
	      perfusion_area_factor, closeNeighborhoodFactor, nu, theta_min,
	      targetICD);
    
    delete gam;    
    delete delta;
    delete eta;
        
    return 0;
} 
