#include<vtkSmartPointer.h>
#include<vtkImageData.h>
#include<vtkPolyData.h>
#include<vtkXMLPolyDataReader.h>
#include<vtkBox.h>
#include<vtkClipPolyData.h>
#include<vtkXMLImageDataWriter.h>
#include<vtkCell.h>
#include<vtkNew.h>

void CoordinatesToPixel(double *point, int *dimensions, int* pixel)
{
  double dx {0.3/dimensions[0]}, dy {0.3/dimensions[1]};
  pixel[0] = floor(dimensions[0]/2.0) + int(point[0]/dx);
  pixel[1] = floor(dimensions[1]/2.0) + int(point[1]/dy);
  pixel[2] = 0;  
}

int main()
{
  int n = 800;			// Number of pixels per axis
  double xmin {-0.15}, xmax {0.15}; // Field of view
  double ymin {-0.15}, ymax {0.15}; // Field of view
  
  
  // Read the source file.
  vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName("..//perfused_geometries/Baseline.vtp");
  reader->Update();  // Needed because of GetScalarRange
  vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();

  // Clipping
  vtkSmartPointer<vtkBox> boundingBox = vtkSmartPointer<vtkBox>::New();
  boundingBox->SetBounds(xmin, xmax,
			 ymin, ymax,
			 -0.15, 0.15);

  vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
  clipper->SetInputData(polyData);
  clipper->SetClipFunction(boundingBox);
  clipper->SetInsideOut(true);
  clipper->Update();
  polyData = clipper->GetOutput();

  vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();
  double *bounds = polyData->GetBounds();
  int dimensions[3] = {n, n, 1};
  int *dims = new int(3);
  
  whiteImage->SetDimensions(dimensions);
  whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR, 1);  

  dims = whiteImage->GetDimensions();
  int z = 0;
  unsigned char *pixels = static_cast<unsigned char *>(whiteImage->GetScalarPointer());
  int extent[6];
  whiteImage->GetExtent(extent);
  cout << "Image extent: ";
  for (const auto& e : extent)
    cout << e << ' ';
  cout << endl;
  cout << "Image dims: ";
  for (int i = 0; i < 3; i++)
    cout << dims[i] << ' ';
  cout << endl;

  
  for(unsigned int x = 0; x < dims[0]; x++)
    {
      for(unsigned int y = 0; y < dims[1]; y++)
  	{
  	  pixels[y * dims[0] + x] = 1;   
  	}
    }

  vtkSmartPointer<vtkPoints> points = polyData->GetPoints();
  vtkSmartPointer<vtkCellArray> lines = polyData->GetLines();
  vtkIdType *indices;
  vtkIdType numberOfPoints;
  unsigned int lineCount = 0;
  for (lines->InitTraversal(); lines->GetNextCell(numberOfPoints, indices); lineCount++)
    {
      // Fill the pixels on the line with 0
      double x0[3], x1[3];
      points->GetPoint(indices[0], x0);
      points->GetPoint(indices[1], x1);
      double dt = 1e-4;
      for (double t = 0; t<=1; t+=dt)
	{
	  double x[3] = {(1-t)*x0[0] + t*x1[0], (1-t)*x0[1] + t*x1[1], 0};
	  int pixel[3];
	  CoordinatesToPixel(x, dims, pixel);
	  pixels[pixel[1]*dims[0] + pixel[0]] = 0;
	}
    }
  
  vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
  writer->SetFileName("Test.vti");
  writer->SetInputData(whiteImage);
  writer->Write();

  return EXIT_SUCCESS;
}
