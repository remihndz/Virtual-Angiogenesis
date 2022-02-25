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

using namespace std;

void Vascularise(string output_filename, string root_tree_filename, string Hull,
		 vector<string> NVR_1, vector<string> NVR_2,
		 long long int n_term_1, long long int n_term_2,
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
  GeneratorData *gen_data = new GeneratorData(160, N_fail, l_lim_fr, perfusion_area_factor,
					      close_neighborhood_factor, 0.25, Delta_nu, 0, false,
					      costEstimator);
					      

  // Domain definition for stage 1
  DomainNVR *domain_1 = new DomainNVR(Hull, NVR_1, n_draw, seed, gen_data);
  domain_1->setIsConvexDomain(true);
  domain_1->setMinBifurcationAngle(theta_min);
  cout << "Domain 1 generated." << endl;

  // Domain definition for stage 2
  DomainNVR *domain_2 = new DomainNVR(Hull, NVR_2, n_draw, seed, gen_data);
  domain_2->setIsConvexDomain(true);
  domain_2->setMinBifurcationAngle(theta_min);
  cout << "Domain 2 generated." << endl;

  StagedDomain *staged_domain = new StagedDomain();
  staged_domain->addStage(n_term_1, domain_1);
  staged_domain->addStage(n_term_2, domain_2);
  cout << "Staged domain initiated." << endl;

  // Checking that the root tree's .cco file exists
  ifstream is_root_tree_correct {root_tree_filename};
  if (!is_root_tree_correct){
    cerr << "Error: file could not be opened. The root tree's .cco file could not be found." << endl;
    exit(1);
  }
  is_root_tree_correct.close();
  
  SingleVesselCCOOTree *tree = new SingleVesselCCOOTree(root_tree_filename, gen_data, gam, delta, eta);
  tree->setIsInCm(true);
  tree->setCurrentStage(2);
  int currentStage{tree->getCurrentStage()};
  cout << "Root tree successfully loaded... ";
  cout << "Current stage is: " << currentStage << endl; 
  cout << "Saving root tree.\n" << endl;
  // Save the root tree
  VTKObjectTreeNodalWriter *tree_writer = new VTKObjectTreeNodalWriter();
  tree_writer->write(output_filename + "_root.vtp", tree);
  tree->save(output_filename + "_root.cco");
  staged_domain->setInitialStage(0);

  long long int n_term_total = n_term_1 + n_term_2 + tree->getNTerms();
  
  StagedFRROTreeGenerator *tree_generator = new StagedFRROTreeGenerator(staged_domain, tree,
									n_term_total,
									{gam, gam},
									{delta, delta},
									{eta, eta});
  tree_generator->setDLim(tree_generator->getDLim()/2.);
  
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
  getline(config, line);
  long long int n_term_1 {stoll(line)};
  getline(config, line);
  long long int n_term_2 {stoll(line)};
  
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
  
  Vascularise(output_filename, input_cco, Hull, NVR_1, NVR_2, n_term_1, n_term_2,
	      gam, delta, eta, n_draw, seed, N_fail, l_lim_fr,
	      perfusion_area_factor, close_neighbourhood_factor, Delta_nu, theta_min);
    
    delete gam;    
    delete delta;
    delete eta;
        
    return 0;
} 
