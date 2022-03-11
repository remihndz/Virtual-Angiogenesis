// STD Libs
#include<cmath>
#include<string>
#include<ctime>
#include<cstdio>
#include<cstdlib>
#include<vector>
#include<fstream>
#include <sstream>

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
#include<structures/domain/NormalDistributionGenerator.h>
#include<structures/vascularElements/SingleVessel.h>
#include<structures/tree/AdimSproutingVolumetricCostEstimator.h>
#include<structures/tree/SproutingVolumetricCostEstimator.h>

using namespace std;

void Vascularise(string output_filename, string root_tree_filename, string Hull,
		 vector<string> NVR_1, vector<string> NVR_2,
		 vector<long long int> nTerms,
		 AbstractConstraintFunction<double, int>* gam,
		 AbstractConstraintFunction<double, int>* delta,
		 AbstractConstraintFunction<double, int>* eta,
		 int n_draw, int seed, int N_fail, double l_lim_fr,
		 double perfusion_area_factor, double close_neighborhood_factor, int Delta_nu,
		 double theta_min)
{
  
  // Simulation parameters for each stage
  SproutingVolumetricCostEstimator *FSprout = new SproutingVolumetricCostEstimator(50, 0.5, 1e+4);
  AbstractCostEstimator *costEstimator = FSprout;
  GeneratorData *gen_data_0 = new GeneratorData(16000, N_fail, l_lim_fr, perfusion_area_factor,
						close_neighborhood_factor, 0.25, Delta_nu, 0, false, costEstimator);

  // Stage for smaller capillaries, centered around the fovea
  GeneratorData *gen_data_1 = new GeneratorData(16000, N_fail, l_lim_fr*0.9, perfusion_area_factor,
						close_neighborhood_factor, 0.25, Delta_nu, 0, false);

					      
  // First stages, large vessels in the parafovea
  StagedDomain *staged_domain = new StagedDomain();
  vector<DomainNVR*> domains;
  for (int i = 0; i<nTerms.size()-1; i++)
    {
      DomainNVR *domain = new DomainNVR(Hull, NVR_1, n_draw, seed, gen_data_0);
      domain->setIsConvexDomain(true);
      domain->setMinBifurcationAngle(theta_min);
      cout << "Domain " << i << " generated." << endl;
      domains.push_back(domain);
    }
  for (int i = 0; i < domains.size(); i++)
    staged_domain->addStage(nTerms[i], domains[i]);
  
  // Domain definition for stage 2
  DomainNVR *fovealDomain = new DomainNVR(Hull, NVR_2, n_draw, seed, gen_data_1);
  fovealDomain->setIsConvexDomain(false);
  fovealDomain->setMinBifurcationAngle(theta_min);
  cout << "Domain " << nTerms.size() << " generated." << endl;
  staged_domain->addStage(nTerms.back(), fovealDomain);

  cout << "Staged domain initialised." << endl;

  // Checking that the root tree's .cco file exists
  ifstream is_root_tree_correct {root_tree_filename};
  if (!is_root_tree_correct){
    cerr << "Error: file could not be opened. The root tree's .cco file could not be found." << endl;
    exit(1);
  }
  is_root_tree_correct.close();
  
  SingleVesselCCOOTree *tree = new SingleVesselCCOOTree(root_tree_filename, gen_data_0, gam, delta, eta);
  tree->setIsInCm(true);
  int currentStage{tree->getCurrentStage()};
  cout << "Current stage is: " << currentStage << endl; 
  cout << "Saving root tree." << endl;
  // Save the root tree
  VTKObjectTreeNodalWriter *tree_writer = new VTKObjectTreeNodalWriter();
  tree_writer->write(output_filename + "_root.vtp", tree);
  tree->save(output_filename + ".root.cco");

  long long int nTermsTotal = tree->getNTerms();
  for (auto& n : nTerms)
    nTermsTotal += n;

  vector<AbstractConstraintFunction<double, int>*> gams, deltas, etas;
  for (int i = 0; i < nTerms.size(); i++)
    {
      gams.push_back(gam);
      deltas.push_back(delta);
      etas.push_back(eta);
    }
  StagedFRROTreeGenerator *tree_generator = new StagedFRROTreeGenerator(staged_domain, tree,
									nTermsTotal,
									gams,
									deltas,
									etas);
  cout << "Staged tree generator initialised." << endl;
  
  cout << "Starting tree generation." << endl;
  tree = {(SingleVesselCCOOTree *) tree_generator->resume(200, "./")};
  cout << "Finished generating the tree." << endl;

  cout << "Saving the results..." << endl;
  tree->save(output_filename + ".cco");  
  tree_writer->write(output_filename + ".vtp", tree);
  cout << "Output written in " << output_filename << ".cco and " << output_filename << ".vtp." << endl;
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
  string output_filename {line};
  
  // Hull
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  string Hull {line};
  cout << "Hull file: " << Hull << endl;				       

  // Domain 1
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  int nb_NVR_1 {stoi(line)};  
  vector<string> NVR_1;
  cout << "The " << nb_NVR_1 << " non vascular regions for domain 1 are: ";
  for (int i = 0; i < nb_NVR_1; i++){
    getline(config, line);
    NVR_1.push_back(line);
    cout << NVR_1[i] << " ";
  }
  cout << endl;

  // Domain 2
  config.ignore(numeric_limits<streamsize>::max(), '\n');
  getline(config, line);
  int nb_NVR_2 {stoi(line)};
  vector<string> NVR_2;
  cout << "The " << nb_NVR_2 << " non vascular regions for domain 2 are: ";
  for (int i = 0; i < nb_NVR_2; i++){
    getline(config, line);
    NVR_2.push_back(line);
    cout << NVR_2[i] << " ";
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
  vector<long long int> nTerms;
  getline(config, line);
  cout << "Reading the number of terminals..." << endl;
  istringstream ss(line);
  long long int n;
  while (ss >> n)		// Read all stages' number of terminal vector (written on one line, separated by blank spaces)
    nTerms.push_back(n);
  cout << "Done." << endl;

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

  // Consecutive attempts to generate a point - N_fail
  int N_fail = 200;
  // Discretisation of the testing triangle for the bifurcation - Delta nu - Figure 1
  int Delta_nu = 7;
  // Buffer size for random point generation
  int n_draw {10000};
  // Random seed
  long long int seed {time(nullptr)};
  
  Vascularise(output_filename, input_cco, Hull, NVR_1, NVR_2, nTerms,
	      gam, delta, eta, n_draw, seed, N_fail, l_lim_fr,
	      perfusion_area_factor, close_neighbourhood_factor, Delta_nu, theta_min);
    
    delete gam;    
    delete delta;
    delete eta;
        
    return 0;
} 
