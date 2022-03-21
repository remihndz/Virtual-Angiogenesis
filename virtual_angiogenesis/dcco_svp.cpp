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
#include<structures/tree/SproutingVolumetricCostEstimator.h>

#include<stats/ObjectTreeStatsManager.h>
#include<stats/VesselObjectHandler.h>

using namespace std;

void Vascularise(string outputFileName, string rootTreeFileName, string Hull,
		 vector<string> NVR1, vector<string> NVR2,
		 long long int nTerms_1, long long int nTerms_2,
		 AbstractConstraintFunction<double, int>* gam,
		 AbstractConstraintFunction<double, int>* delta,
		 AbstractConstraintFunction<double, int>* eta,
		 int nDraw, int seed, int nFail, double l_lim_fr,
		 double perfusion_area_factor, double close_neighborhood_factor, int Delta_nu,
		 double theta_min)
{
  
  // Simulation parameters for each stage
  // SproutingVolumetricCostEstimator *FSprout = new SproutingVolumetricCostEstimator(50, 0.5, 1e+4);
  VolumetricCostEstimator *FSprout = new VolumetricCostEstimator();
  AbstractCostEstimator *costEstimator = FSprout;
  GeneratorData *genData = new GeneratorData(160, nFail, l_lim_fr, perfusion_area_factor,
					      close_neighborhood_factor, 0.25, Delta_nu, 0, false,
					      costEstimator);
					      

  // Domain definition for stage 1
  DomainNVR *domain_1 = new DomainNVR(Hull, NVR1, nDraw, seed, genData);
  domain_1->setIsConvexDomain(true);
  domain_1->setMinBifurcationAngle(theta_min);
  cout << "Domain 1 generated." << endl;

  // Domain definition for stage 2
  DomainNVR *domain_2 = new DomainNVR(Hull, NVR2, nDraw, seed, genData);
  domain_2->setIsConvexDomain(true);
  domain_2->setMinBifurcationAngle(theta_min);
  cout << "Domain 2 generated." << endl;

  StagedDomain *stagedDomain = new StagedDomain();
  stagedDomain->addStage(nTerms_1, domain_1);
  stagedDomain->addStage(nTerms_2, domain_2);
  cout << "Staged domain initialized." << endl;

  // Checking that the root tree's .cco file exists
  ifstream is_root_tree_correct {rootTreeFileName};
  if (!is_root_tree_correct){
    cerr << "Error: file could not be opened. The root tree's .cco file could not be found." << endl;
    exit(1);
  }
  is_root_tree_correct.close();
  
  SingleVesselCCOOTree *tree = new SingleVesselCCOOTree(rootTreeFileName, genData, gam, delta, eta);
  tree->setIsInCm(true);
  tree->setCurrentStage(2);
  int currentStage{tree->getCurrentStage()};
  cout << "Root tree successfully loaded... ";
  cout << "Current stage is: " << currentStage << endl; 
  cout << "Saving root tree.\n" << endl;
  // Save the root tree
  VTKObjectTreeNodalWriter *tree_writer = new VTKObjectTreeNodalWriter();
  tree_writer->write(outputFileName + "_root.vtp", tree);
  tree->save(outputFileName + ".cco.root");
  stagedDomain->setInitialStage(0);

  long long int nTermTotal = nTerms_1 + nTerms_2 + tree->getNTerms();
  
  StagedFRROTreeGenerator *treeGenerator = new StagedFRROTreeGenerator(stagedDomain, tree,
									nTermTotal,
									{gam, gam},
									{delta, delta},
									{eta, eta});
  treeGenerator->setDLim(treeGenerator->getDLim()/2.);
  
  cout << "Staged tree generator initialised." << endl;
  
  cout << "Starting tree generation." << endl;
  tree = {(SingleVesselCCOOTree *) treeGenerator->resume(200, "./")};
  cout << "Finished generating the tree." << endl;

  cout << "Saving the results..." << endl;

  tree->save(outputFileName + "_0.cco");  
  tree_writer->write(outputFileName + "_0.vtp", tree);


  vector<vector<double>> metricsObserver;
  {
    vtkSmartPointer<vtkPolyData> treePolyData = tree->getVtkTree();
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
  }
  cout << "ICD at baseline: " << metricsObserver[0][4] << endl;

  // Domain including perifovea and parafovea but not the FAZ
  DomainNVR *domain = new DomainNVR(Hull, {NVR2.back()}, nDraw, seed, genData);  
  domain->setIsConvexDomain(false);
  domain->setMinBifurcationAngle(theta_min);
  // Add additional vessels until criterion is reached
  double targetICD = 0.022;
  int iter = 1;
  while (metricsObserver[metricsObserver.size()-1][4] < targetICD)
    {
      cout << "Current ICD: " << metricsObserver[metricsObserver.size()-1][4] << endl << endl;
      int n = 50;		// Increment to the number of terminal vessels
      delete stagedDomain;
      stagedDomain = new StagedDomain();
      stagedDomain->addStage(n, domain);
      nTermTotal += n;
      delete treeGenerator;
      treeGenerator = new StagedFRROTreeGenerator(stagedDomain, tree, nTermTotal,
						  {gam}, {delta}, {eta});
	
      tree = {(SingleVesselCCOOTree *) treeGenerator->resume(200, "./")};

      string fileName = outputFileName + to_string(iter) + ".cco";
      tree->save(fileName);
      cout << "Saving checkpoint in " << fileName << endl;
      
      vtkSmartPointer<vtkPolyData> treePolyData = tree->getVtkTree();
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
    }

  // Get the diameter by branching order
  ObjectTreeStatsManager *statsManager = new ObjectTreeStatsManager(tree);
  vector<double> levels, means, stds;
  VesselObjectHandler::ATTRIBUTE att = VesselObjectHandler::ATTRIBUTE::DIAMETER;

  // Saves the tree stats data
  statsManager->getMeanPerLevel(&levels, &means, &stds, att);
  string fileName = outputFileName + "_diameters.dat";
  ofstream outputFile(fileName);
  for (int i = 0; i<levels.size(); i++)
    {
      outputFile << levels[i] << ' ' << means[i] << ' ' << stds[i] << endl;
    }
  outputFile.close();

  // Save the metrics at each stage
  fileName = outputFileName + "_metrics.dat";
  outputFile.open(fileName);
  outputFile << "# NTerms DLim TreeCost NVessels ICD VAD VSD VPI VCI VDI" << endl;
  for (auto e: metricsObserver)
    {
      for (int j = 0; j < e.size(); j++)
	outputFile << e[j] << ' ';
      outputFile << endl;
    }
  outputFile.close();

  cout << "Output written in " << outputFileName << "_i.cco and " << outputFileName << "_i.vtp." << endl;
}







int main(int argc, char *argv[])
{
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
  string outputFileName {line};
  
  // Hull
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  string Hull {line};
  cout << "Hull file: " << Hull << endl;				       

  // Domain 1
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  int nb_NVR1 {stoi(line)};  
  vector<string> NVR1;
  cout << "The " << nb_NVR1 << " non vascular regions for domain 1 are: ";
  for (int i = 0; i < nb_NVR1; i++){
    getline(config, line);
    NVR1.push_back(line);
    cout << NVR1[i] << " ";
  }
  cout << endl;

  // Domain 2
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  int nb_NVR2 {stoi(line)};
  vector<string> NVR2;
  cout << "The " << nb_NVR2 << " non vascular regions for domain 2 are: ";
  for (int i = 0; i < nb_NVR2; i++){
    getline(config, line);
    NVR2.push_back(line);
    cout << NVR2[i] << " ";
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
  long long int nTerms_1 {stoll(line)};
  getline(config, line);
  long long int nTerms_2 {stoll(line)};
  
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  double l_lim_fr {stod(line)};

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
  double close_neighbourhood_factor {stod(line)};

  
  // Copying the config file in the output folder
  string root_results {outputFileName.substr(0, outputFileName.find_last_of("/"))};
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

  // Consecutive attempts to generate a point - nFail
  int nFail = 200;
  // Discretisation of the testing triangle for the bifurcation - Delta nu - Figure 1
  int Delta_nu = 7;
  // Buffer size for random point generation
  int nDraw {10000};
  // Random seed
  long long int seed {time(nullptr)};
  
  Vascularise(outputFileName, inputCCO, Hull, NVR1, NVR2, nTerms_1, nTerms_2,
	      gam, delta, eta, nDraw, seed, nFail, l_lim_fr,
	      perfusion_area_factor, close_neighbourhood_factor, Delta_nu, theta_min);
    
    delete gam;    
    delete delta;
    delete eta;
        
    return 0;
} 
