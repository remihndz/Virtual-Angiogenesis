// VTK Libs
#include <vtkMath.h>
#include <vtkPolygon.h>
#include <vtkRegularPolygonSource.h>
#include <vtkCellArray.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataWriter.h>
#include <iostream>
#include <stdlib.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkClipClosedSurface.h>
#include <vtkDataSetMapper.h>
#include <vtkNamedColors.h>
#include <vtkPlane.h>
#include <vtkPlaneCollection.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkDiskSource.h>

// std Libs
#include<cstdlib>


int CreateSquareWithCircularHole(std::string fileName, double side, double radius)
{

  vtkNew<vtkNamedColors> colors;
  // PolyData to process
  vtkSmartPointer<vtkPolyData> polyData;

  vtkNew<vtkDiskSource> disk;
  disk->SetInnerRadius(radius);
  disk->SetOuterRadius(2*side);
  disk->SetCircumferentialResolution(50);
  disk->Update();

  polyData = disk->GetOutput();
  
  auto center = polyData->GetCenter();
  vtkNew<vtkPlane> plane1;
  plane1->SetOrigin(0, -side/2.0, 0);
  plane1->SetNormal(0, 1, 0);
  vtkNew<vtkPlane> plane2;
  plane2->SetOrigin(0, side/2.0, 0);
  plane2->SetNormal(0, -1, 0);
  vtkNew<vtkPlane> plane3;
  plane3->SetOrigin(-side/2.0, 0, 0);
  plane3->SetNormal(1, 0, 0);
  vtkNew<vtkPlane> plane4;
  plane4->SetOrigin(side/2.0, 0, 0);
  plane4->SetNormal(-1, 0, 0);
  
  vtkNew<vtkPlaneCollection> planes;
  planes->AddItem(plane1);
  planes->AddItem(plane2);
  planes->AddItem(plane3);
  planes->AddItem(plane4);

  vtkNew<vtkClipClosedSurface> clipper;
  clipper->SetInputData(polyData);
  clipper->SetClippingPlanes(planes);
  clipper->SetActivePlaneId(2);
  clipper->SetScalarModeToColors();
  clipper->SetClipColor(colors->GetColor3d("Banana").GetData());
  clipper->SetBaseColor(colors->GetColor3d("Tomato").GetData());
  clipper->SetActivePlaneColor(colors->GetColor3d("SandyBrown").GetData());

  vtkNew<vtkDataSetMapper> clipMapper;
  clipMapper->SetInputConnection(clipper->GetOutputPort());

  vtkNew<vtkActor> clipActor;
  clipActor->SetMapper(clipMapper);
  clipActor->GetProperty()->SetColor(1.0000, 0.3882, 0.2784);
  clipActor->GetProperty()->SetInterpolationToFlat();

  // // Create graphics stuff
  vtkNew<vtkRenderer> ren1;
  ren1->SetBackground(colors->GetColor3d("SteelBlue").GetData());

  vtkNew<vtkRenderWindow> renWin;
  renWin->AddRenderer(ren1);
  renWin->SetSize(512, 512);
  renWin->SetWindowName("ClipClosedSurface");

  vtkNew<vtkRenderWindowInteractor> iren;
  iren->SetRenderWindow(renWin);

  // Add the actors to the renderer, set the background and size
  ren1->AddActor(clipActor);

  // Generate an interesting view
  ren1->ResetCamera();
  ren1->GetActiveCamera()->Azimuth(120);
  ren1->GetActiveCamera()->Elevation(30);
  ren1->GetActiveCamera()->Dolly(1.0);
  ren1->ResetCameraClippingRange();

    // Save the polydata
  vtkNew<vtkPolyDataWriter> writer;
  writer->SetFileName(fileName.c_str());
  writer->SetInputData(clipper->GetOutput());
  writer->Write();
  
  return EXIT_SUCCESS;
}

int CreateVascularAndAvascularRegions(double FOV_in_cm, double FAZ_in_cm)
{

  // Create the hull, i.e., the square representing the OCTA FOV 
  {
    vtkNew<vtkPoints> points;
    points->InsertNextPoint(-FOV_in_cm/2., -FOV_in_cm/2., 0.0);
    points->InsertNextPoint(FOV_in_cm/2., -FOV_in_cm/2., 0.0);
    points->InsertNextPoint(FOV_in_cm/2., FOV_in_cm/2., 0.0);
    points->InsertNextPoint(-FOV_in_cm/2., FOV_in_cm/2., 0.0);

    std::cout << "Points created." << std::endl;
  
    // Polygon and CellArray to store it

    vtkNew<vtkCellArray> aCellArray;
    vtkNew<vtkPolygon> aPolygon;
  
    aPolygon->GetPointIds()->SetNumberOfIds(4);
    aPolygon->GetPointIds()->SetId(0,0);
    aPolygon->GetPointIds()->SetId(1, 1);
    aPolygon->GetPointIds()->SetId(2, 2);
    aPolygon->GetPointIds()->SetId(3, 3);

    std::cout << "Polygon generated." << std::endl;
  
    aCellArray->InsertNextCell(aPolygon);

    std::cout << "Cell Array filled." << std::endl;
  
    // PolyData to write result
    vtkNew<vtkPolyData> aPolyData;
    aPolyData->SetPoints(points);
    aPolyData->SetPolys(aCellArray);
  
    std::cout << "Polydata created, moving on to writing it." << std::endl;
  
    // Save the polydata
    vtkNew<vtkPolyDataWriter> writer;
    std::string fileName = "../../base_geometries/VTKFiles/hull.vtk";
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(aPolyData);
    writer->Write();
    std::cout << "Hull created and saved in " << fileName << std::endl;
  }
      
  // Create NonVascular Regions for stage 1
  {
    vtkNew<vtkRegularPolygonSource> disk;
    disk->SetNumberOfSides(50);
    disk->SetRadius(FOV_in_cm/3.0);
    disk->Update();  
    std::string fileName = "../../base_geometries/VTKFiles/non_vascular_region_1.vtk";

    vtkNew<vtkPolyDataWriter> writer;
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(disk->GetOutput());
    writer->Write();
    std::cout << "Non Vascular Region for stage 1 written in " << fileName << std::endl;
  }
  
  // Create NonVascular Regions for stage 2
  {
    vtkNew<vtkRegularPolygonSource> disk;
    disk->SetNumberOfSides(50);
    disk->SetRadius(FAZ_in_cm);
    disk->SetCenter(0,0,0);
    disk->Update();  

    std::string fileName = "../../base_geometries/VTKFiles/non_vascular_region_2_inner.vtk";
    vtkNew<vtkPolyDataWriter> writer;
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(disk->GetOutput());
    writer->Write();
    std::cout << "Inner Non Vascular Region for stage 2 written in " << fileName << std::endl;

    fileName = "../../base_geometries/VTKFiles/non_vascular_region_2_outer.vtk";
    CreateSquareWithCircularHole(fileName, FOV_in_cm, FOV_in_cm/3.0);
    std::cout << "Outer Non Vascular Region for stage 2 written in " << fileName << std::endl;
  }
}

int main()
{
  std::cout << "Starting creation of vascular and non vascular regions for a 2 stages CCO growth..." << std::endl;

  const int dir_err = std::system("mkdir -p ../../base_geometries/VTKFiles/");
  if (-1 == dir_err)
    {
      printf("Error creating directory!n");
      exit(1);
    }

  double FOV_in_cm = 0.3;	// Field of View, typically 6mm or 3mm
  double FAZ_in_cm = 0.05;	// Foveal Avascular Zone, typically the FAZ has a diameter of 0.5mm in the SVP
  
  CreateVascularAndAvascularRegions(FOV_in_cm, FAZ_in_cm);

  return EXIT_SUCCESS;
}

    
    
