#include<cmath>

#include<vtkSmartPointer.h>
#include<vtkImageData.h>
#include<vtkPolyData.h>
#include<vtkXMLPolyDataReader.h>
#include<vtkBox.h>
#include<vtkClipPolyData.h>
#include<vtkXMLImageDataWriter.h>
#include<vtkCell.h>
#include<vtkNew.h>
#include<vtkImageEuclideanDistance.h>
#include<vtkImageShiftScale.h>
#include <vtkImageDilateErode3D.h>

void CoordinatesToPixel(double *point, int *dimensions, int* pixel)
{
  double dx {0.3/dimensions[0]}, dy {0.3/dimensions[1]};
  pixel[0] = floor(dimensions[0]/2.0) + int(point[0]/dx);
  pixel[1] = floor(dimensions[1]/2.0) + int(point[1]/dy);
  pixel[2] = 0;  
}

void PixelToCoordinates(int *pixel, int *dimensions, double* point)
{
  double dx {0.3/dimensions[0]}, dy {0.3/dimensions[1]};
  point[0] = -0.15 + dx * pixel[0];
  point[1] = -0.15 + dy * pixel[1];
}

int main()
{
  int P;
  cout << "Input the power P for the number of pixels (nbPixel = 2^P)" << endl;
  cin >> P;
  int n = pow(2, P);			// Number of pixels per axis
  double xmin {-0.15}, xmax {0.15}; // Field of view
  double ymin {-0.15}, ymax {0.15}; // Field of view
    
  // Read the source file.
  vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName("../Results/Results50Sims/sim1.vtp");
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

  // Print the polydata tree on a vtkImageData where background is 1 and tree is 0
  vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
  double *bounds = polyData->GetBounds();
  int dimensions[3] = {n, n, 1};
  int *dims = new int(3);
  double spacing[3] = {(bounds[1]-bounds[0])/n, (bounds[3]-bounds[2])/n, 1};
  
  imageData->SetDimensions(dimensions);
  imageData->AllocateScalars(VTK_UNSIGNED_CHAR, 1);  

  dims = imageData->GetDimensions();
  int z = 0;
  unsigned char *pixels = static_cast<unsigned char *>(imageData->GetScalarPointer());
  int extent[6];
  imageData->GetExtent(extent);
  cout << "Image extent: ";
  for (const auto& e : extent)
    cout << e << ' ';
  cout << endl;
  cout << "Image dims: ";
  for (int i = 0; i < 3; i++)
    cout << dims[i] << ' ';
  cout << endl;

  
  for(int i = 0; i < dims[0]; i++)
    {
      for(int j = 0; j < dims[1]; j++)
  	{
	  double point[3];
	  int pixel[3] = {i,j,0};
	  PixelToCoordinates(pixel, dims, point);
	  if (pow(point[0],2) + pow(point[1],2) < pow(0.04 ,2))
	    pixels[i * dims[0] + j] = 0;
	  else
	    pixels[i * dims[0] + j] = 1;   
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
      double dt = 1e-5;
      for (double t = 0; t<=1; t+=dt)
	{
	  double x[3] = {(1-t)*x0[0] + t*x1[0], (1-t)*x0[1] + t*x1[1], 0};
	  int pixel[3];
	  CoordinatesToPixel(x, dims, pixel);
	  pixels[pixel[1]*dims[0] + pixel[0]] = 0;
	}
    }
  
  vtkSmartPointer<vtkImageDilateErode3D> dilateErode = vtkSmartPointer<vtkImageDilateErode3D>::New();
  dilateErode->SetInputData(imageData);
  dilateErode->SetDilateValue(0);
  dilateErode->SetErodeValue(1);
  dilateErode->SetKernelSize(2,2,1);
  dilateErode->Update();
  imageData = dilateErode->GetOutput();

  double *spacingDT = imageData->GetSpacing();
  cout << "Dilate Filter spacings: " << spacingDT[0] << ' ' << spacingDT[1] << ' ' << spacingDT[2] << endl;  
  vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
  writer->SetFileName("Test.vti");
  writer->SetInputData(imageData);
  writer->Write();


  // Compute the euclidean distance map
  cout << "Spacings: " << spacing[0] << ' ' << spacing[1] << ' ' << spacing[2] << endl;  
  vtkSmartPointer<vtkImageEuclideanDistance> dtFilter = vtkSmartPointer<vtkImageEuclideanDistance>::New();
  dtFilter->SetAlgorithmToSaitoCached();
  dtFilter->SetDimensionality(2);
  dtFilter->SetInputData(imageData);
  dtFilter->Update();

  spacingDT = dtFilter->GetOutput()->GetSpacing();
  cout << "DT Filter spacings: " << spacingDT[0] << ' ' << spacingDT[1] << ' ' << spacingDT[2] << endl;

  cout << "Created the filter" << endl;
  
  vtkSmartPointer<vtkImageData> distanceMap = vtkSmartPointer<vtkImageData>::New();

  cout << "Created the distance map" << endl;
  // Scale by the size of a pixel (distanceMap is the distance in pixels at this point)
  vtkSmartPointer<vtkImageShiftScale> scaler = vtkSmartPointer<vtkImageShiftScale>::New();
  scaler->SetOutputScalarTypeToDouble();
  scaler->SetInputConnection(dtFilter->GetOutputPort());
  scaler->SetScale(spacing[0]);
  scaler->Update();

  distanceMap = scaler->GetOutput();
  distanceMap->SetSpacing(spacing);
  spacingDT = distanceMap->GetSpacing();
  cout << "Distance map spacings: " << spacingDT[0] << ' ' << spacingDT[1] << ' ' << spacingDT[2] << endl;

  // Compute average intercapillary distance
  // pixels = static_cast<unsigned char *>(distanceMap->GetScalarPointer());
  double ICD = 0;
  int count = 0;
  for(int i = 0; i < dims[0]; i++)
    {
      for(int j = 0; j < dims[1]; j++)
  	{
	  float dist = distanceMap->GetScalarComponentAsFloat(i,j,0,0);
	  if (dist > 0)
	    {
	      ICD += dist;
	      count++;
	    }
  	}
    }
  
  cout << "ICD: " << ICD/count << " with " << count << " pixels counted." << endl;
  
  writer->SetFileName("DistanceMap.vti");
  writer->SetInputData(distanceMap);
  writer->Write();
  
  
  return EXIT_SUCCESS;
}
