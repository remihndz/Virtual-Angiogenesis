#include<VascularIndexes.hpp>



double VesselSkeletonDensity(const vector<SingleVessel *> Vessels, double domainAreaInCm, int ignoreUpToStage)
{
  double totalSkeletonLength {0.0};
  for (auto v: Vessels)
    {
      if (v->stage >= ignoreUpToStage)
	totalSkeletonLength += v->length; // In cm
    }

  return totalSkeletonLength/domainAreaInCm * 1e-4; // Converting in micron^-1
}

double VesselComplexityIndex(const vector<SingleVessel *>Vessels, double domainAreaInCm, int ignoreUpToStage)
{
  double totalPerimeter {0.0}, totalArea{0.0};
  for (auto v: Vessels)
    {
      if (v->stage >= ignoreUpToStage)
	{
	  totalPerimeter += v->length*2.0; // In cm
	  totalArea      += v->length*2.0*v->radius; // In cm^2
	}
    }

  return pow(totalPerimeter,2)/(4.0*totalArea*M_PI)/1.5711; // Normalized as in Chu 2016
}

double VesselAreaDensity(const vector<SingleVessel *>Vessels, double domainAreaInCm, int ignoreUpToStage)
{
  double totalArea {0.0};
  for (auto v: Vessels)
    {
      if (v->stage >= ignoreUpToStage)
	totalArea += v->length*2.0*v->radius; // In cm^2
    }
  return totalArea/domainAreaInCm;
}

double VesselPerimeterIndex(const vector<SingleVessel *>Vessels, double domainAreaInCm, int ignoreUpToStage)
// Note: this is twice the vessel skeleton density
{
  double totalPerimeter {0.0};
  for (auto v: Vessels)
    {
      if (v->stage >= ignoreUpToStage)
	totalPerimeter += 2*v->length; // In cm
    }
  return totalPerimeter/domainAreaInCm * 1e-4; // In micron^-1
}


double VesselDiameterIndex(const vector<SingleVessel *> Vessels, double domainAreaInCm, int ignoreUpToStage)
{
  double totalLength {0.0}, totalArea {0.0};
  for (auto v: Vessels)
    {
      if (v->stage >= ignoreUpToStage)
	{
	  totalLength += v->length; // In cm
	  totalArea   += 2.0*v->radius*v->length; // In cm^2
	}
    }
  return totalArea/totalLength * 1e4; // In micron
}

