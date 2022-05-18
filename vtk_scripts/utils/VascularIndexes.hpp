/* 
   Computes the vascular indexes proposed by Chu et al. 
   The argument ignoreUpToStage is used to only include
   vessel which were grown at stage ignoreUpToStage or 
   later. Set to -1 to ignore the linking vessels 
   lying outside the perfusion domain.
*/


#include<vector>
#include<ctime>
#include<cmath>
#include<unordered_map>

// VItA Libs
#include<structures/vascularElements/AbstractVascularElement.h>
#include<structures/vascularElements/SingleVessel.h>

double VesselSkeletonDensity(const vector<SingleVessel *> Vessels, double domainAreaInCm=0.09, int ignoreUpToStage = -1);

double VesselComplexityIndex(const vector<SingleVessel *> Vessels, double domainAreaInCm = 0.09, int ignoreUpToStage = -1);

double VesselAreaDensity(const vector<SingleVessel *> Vessels, double domainAreaInCm = 0.09, int ignoreUpToStage = -1);

double VesselPerimeterIndex(const vector<SingleVessel *> Vessels, double domainAreaInCm = 0.09, int ignoreUpToStage = -1);

double VesselDiameterIndex(const vector<SingleVessel *> Vessels, double domainAreaInCm = 0.09, int ignoreUpToStage = -1);
