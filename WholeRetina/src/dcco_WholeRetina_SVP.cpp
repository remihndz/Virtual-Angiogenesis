// STD Libs
#include<cmath>
#include<string>
#include<ctime>
#include<cstdio>
#include<cstdlib>
#include<vector>
#include<fstream>
#include<ctime>
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
#include<structures/domain/SimpleDomain.h>
#include<structures/domain/DummyDomain.h>
#include<structures/domain/NormalDistributionGenerator.h>
#include<structures/vascularElements/SingleVessel.h>
#include<structures/tree/AdimSproutingVolumetricCostEstimator.h>
#include<structures/tree/SproutingVolumetricCostEstimator.h>

using namespace std;


string resultsDir = "../Results/";
string outputFileName = resultsDir + "110terminals",
  rootVesselCCO = "rootVessel.cco";
long long int n_term_0 {30}, n_term_1 {100}; // Number of terminals to be generated in stage 0 and 1

void vascularise(string output_filename, AbstractConstraintFunction<double, int>* gam,
    AbstractConstraintFunction<double, int>* delta, AbstractConstraintFunction<double, int>* eta,
    int n_draw, int seed, int N_fail, double l_lim_fr,
    double perfusion_area_factor, double close_neighborhood_factor, int Delta_nu,
    double theta_min)
{

    // Boundary conditions - Flow, radius and position of the inlet 
    double q0 {15.0};		// muL/s
    double r0 {0.0070};		// cm
    point x0 {0.0, 0.0, -1.0};	// cm
    double p0 {25 * 133.3224};  // mmHg*133.3224 = Pa, should be consistent with pressure unit 
    
    // Domain geometry
    string vtkRetina = "../vtkFiles/Results/WholeRetina.vtk";
    string vtkFAZ    = "../vtkFiles/Results/FAZ.vtk";

    string vtkRetinaHull = "../vtkFiles/Results/RetinaHull.vtk",
      vtkRetinaNVR = "../vtkFiles/Results/RetinaNVR.vtk";
    
    // Optimisation parameters
    GeneratorData *gen_data_0 = new GeneratorData(16000, N_fail, l_lim_fr, perfusion_area_factor, close_neighborhood_factor, 0.1, Delta_nu, 1, false);
    GeneratorData *gen_data_1 = new GeneratorData(16000, N_fail, l_lim_fr, perfusion_area_factor, close_neighborhood_factor, 0.1, Delta_nu, 0, false);

    
    // Domain definition
    // SimpleDomain *domain_0 = new SimpleDomain(vtkRetina, n_draw, seed, gen_data_0);
    DomainNVR *domain_0 = new DomainNVR(vtkRetinaHull, {vtkRetinaNVR, vtkFAZ}, n_draw, seed, gen_data_0);
    long long int seed_1 {time(nullptr)};
    DomainNVR *domain_1 = new DomainNVR(vtkRetinaHull, {vtkRetinaNVR, vtkFAZ}, n_draw, seed_1, gen_data_1);
	
    domain_0->setMinBifurcationAngle(theta_min);
    domain_0->setIsBifPlaneContrained(true);
    domain_0->setMinPlaneAngle(theta_min);
    domain_0->setIsConvexDomain(true);
    domain_1->setMinBifurcationAngle(theta_min);
    domain_1->setIsBifPlaneContrained(true);
    domain_1->setMinPlaneAngle(theta_min);
    domain_1->setIsConvexDomain(true);

    SingleVesselCCOOTree *tree = new SingleVesselCCOOTree(x0, r0, q0, gam, delta, eta, p0, 1.0e-5,
							  gen_data_0);    
    VTKObjectTreeNodalWriter *tree_writer = new VTKObjectTreeNodalWriter();
    // SingleVesselCCOOTree *tree = new SingleVesselCCOOTree(rootVesselCCO, gen_data_0, gam, delta, eta);
    // tree->save("../Results/Root.cco");
    // tree_writer->write("../Results/Root.vtp", tree);
    

    // Define domain stages
    tree->setCurrentStage(0);
    StagedDomain *staged_domain = new StagedDomain();
    staged_domain->addStage(n_term_0, domain_0);
    staged_domain->addStage(n_term_1, domain_1);
    
    // Creation of DCCO generator
    StagedFRROTreeGenerator *tree_generator = new StagedFRROTreeGenerator(staged_domain, tree,
									  n_term_0+n_term_1,
									  {gam, gam}, {delta, {new ConstantConstraintFunction<double, int>(0.3)}}, {eta, eta});

    // Indicates that all meshes and boundary conditions are in cm (important for the viscosity computation).
    tree->setIsInCm(true);

    cout << "Finished preparing. Moving on to generation." << endl;
    
    // Executes DCCO generator
    tree = {(SingleVesselCCOOTree *) tree_generator->generate(10, ".")};
    
    // Saves the outputs as CCO and VTP files
    tree->save(output_filename + ".cco");
    tree_writer->write(output_filename + ".vtp", tree);

    cout << "Output written in " << output_filename << ".cco and " << output_filename << ".vtp." << endl;
          
    delete tree_writer;
    delete tree_generator;
    delete staged_domain;
    delete domain_0;
    delete gen_data_0;
}

int main(int argc, char *argv[])
{
    // Consecutive attempts to generate a point - N_fail
    int N_fail = 20;
    // Correction step factor - fr - Eq (14)
    double l_lim_fr = 0.9;
    // Discretisation of the testing triangle for the bifurcation - Delta nu - Figure 1
    int Delta_nu = 8;

    // Geometrical constrains - 
    // Power-law coefficient - gamma - Eq (3) - Murray's law
    AbstractConstraintFunction<double,int> *gam {new ConstantConstraintFunction<double, int>(3.)};
    // Symmetry ratio parameter - delta - Eq (4)
    AbstractConstraintFunction<double,int> *delta {new ConstantConstraintFunction<double, int>(0.8)};
    // Viscosity in cP - eta - Eq (7)
    AbstractConstraintFunction<double,int> *eta {new ConstantConstraintFunction<double, int>(3.6)};

    // Buffer size for random point generation
    int n_draw {10000};
    // Random seed
    long long int seed {time(nullptr)};
    // Minimum bifurcation angle - theta_min - Eq (17) 
    double theta_min {(3./18.) * M_PI}; // {20.0*M_PI/180}
    // l_min tuning parameter - nu - Eq (12)
    double perfusion_area_factor {.9};
    // Heuristic parameter to reduce the neighbour vessels to be tested for connections in line 15 of Algorithm 2
    double close_neighbourhood_factor {4.0};

    vascularise(outputFileName, gam, delta, eta, n_draw, seed, N_fail, l_lim_fr, perfusion_area_factor, close_neighbourhood_factor, Delta_nu, theta_min);
    
    delete gam;    
    delete delta;
    delete eta;
        
    return 0;
} 
