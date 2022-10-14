// VTK Libs
#include <vtkMath.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkNew.h>

// #include <vtkPolygon.h>
// #include <vtkRegularPolygonSource.h>
// #include <vtkCellArray.h>
// #include <vtkPoints.h>
// #include <iostream>
// #include <stdlib.h>
// #include <vtkActor.h>
// #include <vtkCamera.h>
// #include <vtkClipClosedSurface.h>
// #include <vtkDataSetMapper.h>
// #include <vtkNamedColors.h>
// #include <vtkPlane.h>
// #include <vtkPlaneCollection.h>
// #include <vtkProperty.h>
// #include <vtkRenderWindow.h>
// #include <vtkRenderWindowInteractor.h>
// #include <vtkRenderer.h>
// #include <vtkDiskSource.h>

// std Libs
#include<cstdlib>
#include<cmath>

double AreaSphericalCap(double radius, double angle){
  return 2*M_PI*radius*radius*(1-cos(angle));
}

double AngleFromArea(double radius, double targetArea){
  return acos( 1 - targetArea/(2*M_PI*radius*radius) );
}

int CreateRetina(char* fileNameRetina, double radius, double cropAngle,
		 char* fileNameFAZ, double FAZCenter, double FAZArea)
{

  // Creating the polydata of the retina
  
  vtkNew<vtkSphereSource> sphere;
  vtkNew<vtkPolyDataWriter> writer;

  sphere->SetCenter(0.0, 0.0, 0.0);
  sphere->SetRadius(radius);
  sphere->SetPhiResolution(500);
  sphere->SetThetaResolution(500);
  sphere->SetStartPhi(cropAngle/2.0);
  sphere->SetEndPhi(360-cropAngle/2.0);
  sphere->Update();
  
  writer->SetInputData(sphere->GetOutput());
  writer->SetFileName(fileNameRetina);
  // writer->SetFileVersion(42);	// Maybe only necessary in newest versions that write in vtk format 5.1 by default (e.g., in python)
  writer->Write();


  // Polydata of the FAZ, centered at FAZCenter degrees (kappa angle) from the pupillary axis
  // and an area of FAZArea
  double widthAngle = AngleFromArea(radius, FAZArea);
  sphere->SetStartPhi(180+FAZCenter-widthAngle);
  sphere->SetEndPhi(180+FAZCenter+widthAngle);
  sphere->Update();

  writer->SetFileName(fileNameFAZ);
  writer->Write();

  std::cout << "Target FAZ area is " << FAZArea << " and generated FAZ has area " << AreaSphericalCap(radius, widthAngle) << std::endl;
  std::cout << "This yields a disc-fovea distance (considering the disc is on the pupillary axis) of " << radius*FAZCenter*M_PI/180 << std::endl;
  
  return EXIT_SUCCESS;
}


// int CreateFAZ(std::string fileName, double radiusFAZ, double radiusRetina)
// {
//   vtkSmartPointer<vtkSphereSource> sphere;


  
int main()
{
  const int dir_err = std::system("mkdir -p ../Results/");
  if (-1 == dir_err)
    {
      printf("Error creating directory!n");
      exit(1);
    }
  
  double FAZ_in_cm = 0.04;	// Foveal Avascular Zone, typically the FAZ has a diameter of 0.5mm in the SVP
  double radiusRetina = 1.0; 	// Typical radius of a human retina (from Hutton-Smith 2016)
  double cropAngle = 120;
  
  CreateRetina("../Results/WholeRetina.vtk", radiusRetina, cropAngle,
	       "../Results/FAZ.vtk", 20, 0.003);

  return EXIT_SUCCESS;
}

    
    
