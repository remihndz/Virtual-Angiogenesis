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

int CreateVascularAndAvascularRegions()
{
  // Create the hull
  {
    vtkNew<vtkPoints> points;
    points->InsertNextPoint(-3.0, -3.0, 0.0);
    points->InsertNextPoint(3.0, -3.0, 0.0);
    points->InsertNextPoint(3.0, 3.0, 0.0);
    points->InsertNextPoint(-3.0, 3.0, 0.0);

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
    std::string fileName = "../../base_geometries/plane.vtk";
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(aPolyData);
    writer->Write();
    std::cout << "Hull created and saved in " << fileName << std::endl;
  }
    
  // Create Vascular Region for stage 1
  {
    std::string fileName = "../../base_geometries/vascular_region_1.vtk";
    CreateSquareWithCircularHole(fileName, 6.0, 1.);
    std::cout << "Vascular region for stage 1 written in " << fileName << std::endl;
  }
  
  // Create NonVascular Regions for stage 1
  {
    vtkNew<vtkRegularPolygonSource> disk;
    disk->SetNumberOfSides(50);
    disk->SetRadius(2.0);
    disk->Update();  
    std::string fileName = "../../base_geometries/non_vascular_region_1.vtk";

    vtkNew<vtkPolyDataWriter> writer;
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(disk->GetOutput());
    writer->Write();
    std::cout << "Non Vascular Region for stage 1 written in " << fileName << std::endl;
  }

  // Create Vascular Region for stage 2
  {
    vtkNew<vtkDiskSource> disk;
    disk->SetCircumferentialResolution(50);
    disk->SetOuterRadius(1.);
    disk->SetInnerRadius(0.4);
    //disk->SetCenter(0.,0.,0.);
    disk->Update();  
    std::string fileName = "../../base_geometries/vascular_region_2.vtk";

    vtkNew<vtkPolyDataWriter> writer;
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(disk->GetOutput());
    writer->Write();
    std::cout << "Vascular Region for stage 2 written in " << fileName << std::endl; 
  }
  
  // Create NonVascular Regions for stage 2
  {
    vtkNew<vtkRegularPolygonSource> disk;
    disk->SetNumberOfSides(50);
    disk->SetRadius(0.4);
    disk->SetCenter(0,0,0);
    disk->Update();  

    std::string fileName = "../../base_geometries/non_vascular_region_2_inner.vtk";
    vtkNew<vtkPolyDataWriter> writer;
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(disk->GetOutput());
    writer->Write();
    std::cout << "Inner Non Vascular Region for stage 2 written in " << fileName << std::endl;

    fileName = "../../base_geometries/non_vascular_region_2_outer.vtk";
    CreateSquareWithCircularHole(fileName, 6.0, 2.0);
    std::cout << "Outer Non Vascular Region for stage 2 written in " << fileName << std::endl;
  }
}

int main()
{
  std::cout << "Starting creation of vascular and non vascular regions for a 2 stages CCO growth..." << std::endl;

  CreateVascularAndAvascularRegions();

  return EXIT_SUCCESS;
}

    
    
// int main(int argc, char **argv){

//   std::cout << "Starting generation of the plane..." << std::endl;
//   // Vertex of the square
//   vtkNew<vtkPoints> points;
//   points->InsertNextPoint(-3.0, -3.0, 0.0);
//   points->InsertNextPoint(3.0, -3.0, 0.0);
//   points->InsertNextPoint(3.0, 3.0, 0.0);
//   points->InsertNextPoint(-3.0, 3.0, 0.0);

//   std::cout << "Points created." << std::endl;
  
//   // Polygon and CellArray to store it

//   vtkNew<vtkCellArray> aCellArray;
//   vtkNew<vtkPolygon> aPolygon;
  
//   aPolygon->GetPointIds()->SetNumberOfIds(4);
//   aPolygon->GetPointIds()->SetId(0,0);
//   aPolygon->GetPointIds()->SetId(1, 1);
//   aPolygon->GetPointIds()->SetId(2, 2);
//   aPolygon->GetPointIds()->SetId(3, 3);

//   std::cout << "Polygon generated." << std::endl;
  
//   aCellArray->InsertNextCell(aPolygon);

//   std::cout << "Cell Array filled." << std::endl;
  
//   // PolyData to write result
//   vtkNew<vtkPolyData> aPolyData;
//   aPolyData->SetPoints(points);
//   aPolyData->SetPolys(aCellArray);
  
//   std::cout << "Polydata created, moving on to writing it." << std::endl;
  
//   // Save the polydata
//   vtkNew<vtkPolyDataWriter> writer;
//   writer->SetFileName("../../base_geometries/plane.vtk");
//   writer->SetInputData(aPolyData);
//   writer->Write();

//   // Creates three disk corresponding to nonvascularized areas at each stage of the CCO
//   vtkNew<vtkRegularPolygonSource> disk;
//   disk->SetNumberOfSides(50);
//   disk->SetRadius(2.0);
//   disk->SetCenter(0,0,0);
//   disk->Update();  
//   writer->SetFileName("../../base_geometries/NVR_1.vtk");
//   writer->SetInputData(disk->GetOutput());
//   writer->Write();

//   disk->SetRadius(1.0);
//   disk->Update();
//   writer->SetFileName("../../base_geometries/NVR_2.vtk");
//   writer->SetInputData(disk->GetOutput());
//   writer->Write();
  
//   disk->SetRadius(0.4);
//   disk->Update();
//   writer->SetFileName("../../base_geometries/NVR_3.vtk");
//   writer->SetInputData(disk->GetOutput());
//   writer->Write();
  
//   return EXIT_SUCCESS;
// }
