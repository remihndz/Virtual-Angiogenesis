#include<cmath>

// VTK libs
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

// OpenCV libs
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>


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

cv::Mat convertVtkImageDataToCVMat(const vtkSmartPointer<vtkImageData> &vtkImage) {
  int imageDimensions[3] = {0, 0, 0}; // Width, Hight, Depth --> Depth is not equal to number of image channels!
  vtkImage->GetDimensions(imageDimensions);
  int imageWidth = imageDimensions[0];
  int imageHeight = imageDimensions[1];
  int numberOfImageChannels = vtkImage->GetNumberOfScalarComponents();
  int cvType = 0;
  switch(numberOfImageChannels){
    case 1: cvType = CV_8UC1; break;
    case 3: cvType = CV_8UC3; break;
    case 4: cvType = CV_8UC4; break;
    default: std::cerr << "Check number of vtk image channels!" << std::endl;
  }
  auto resultingCVMat = cv::Mat(imageHeight, imageWidth, cvType);
  // Loop over the vtkImageData contents.
  for ( int heightPos = 0; heightPos < imageHeight; heightPos++ ){
    for ( int widthPos = 0; widthPos < imageWidth; widthPos++ ){
      auto pixel = static_cast<unsigned char *>(vtkImage->GetScalarPointer(widthPos, heightPos, 0));
      resultingCVMat.at<unsigned char>(heightPos, widthPos) = *pixel;
    }
  }
  return resultingCVMat;
}


int main()
{
  int P;
  cout << "Input the power P for the number of pixels (nbPixel = 2^P)" << endl;
  cin >> P;
  int n = pow(2, P);			// Number of pixels per axis
  double xmin {-0.15}, xmax {0.15}; // Field of view
  double ymin {-0.15}, ymax {0.15}; // Field of view

  unsigned int foregroundPixel {0}, backgroundPixel {255};
  
  // Read the source file.
  vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName("../Results/Results50Sims/sim1.vtp");
  reader->Update();  // Creates a deprecated warning
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

  // Print the polydata tree on a vtkImageData 
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
	  // if (pow(point[0],2) + pow(point[1],2) < pow(0.04 ,2))
	  //   pixels[i * dims[0] + j] = foregroundPixel;
	  // else
	  pixels[i * dims[0] + j] = backgroundPixel;   
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
	  pixels[pixel[1]*dims[0] + pixel[0]] = foregroundPixel;
	}
    }
  
  vtkSmartPointer<vtkImageDilateErode3D> dilateErode = vtkSmartPointer<vtkImageDilateErode3D>::New();
  dilateErode->SetInputData(imageData);
  dilateErode->SetDilateValue(0);
  dilateErode->SetErodeValue(1);
  dilateErode->SetKernelSize(2,2,1);
  dilateErode->Update();
  imageData = dilateErode->GetOutput();

  vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
  writer->SetFileName("Test.vti");
  writer->SetInputData(imageData);
  writer->Write();


  // Compute the euclidean distance map
  vtkSmartPointer<vtkImageEuclideanDistance> dtFilter = vtkSmartPointer<vtkImageEuclideanDistance>::New();
  dtFilter->SetAlgorithmToSaitoCached();
  dtFilter->SetDimensionality(2);
  dtFilter->SetInputData(imageData);
  dtFilter->Update();
  vtkSmartPointer<vtkImageData> distanceMap = vtkSmartPointer<vtkImageData>::New();

  // Scale by the size of a pixel (distanceMap is the distance in pixels at this point)
  vtkSmartPointer<vtkImageShiftScale> scaler = vtkSmartPointer<vtkImageShiftScale>::New();
  scaler->SetOutputScalarTypeToDouble();
  scaler->SetInputConnection(dtFilter->GetOutputPort());
  scaler->SetScale(spacing[0]);
  scaler->Update();

  distanceMap = scaler->GetOutput();
  distanceMap->SetSpacing(spacing);

  cv::Mat bwImage = convertVtkImageDataToCVMat(imageData);
  // cv::imshow("Black Background Image", bwImage);
  cv::Mat dist;
  cv::distanceTransform(bwImage, dist, cv::DIST_L2, 3); // The last number is the mask size. Use 5 for more precise calculations
  // Normalize the distance image for range = {0.0, 1.0}
  // so we can visualize and threshold it
  cv::normalize(dist, dist, 0, 1.0, cv::NORM_MINMAX);
  // cv::imshow("Distance Transform Image", dist);
  // Threshold to obtain the peaks
  // This will be the markers for the foreground objects
  cv::threshold(dist, dist, 0.4, 1.0, cv::THRESH_BINARY);
  // Dilate a bit the dist image
  cv::Mat kernel1 = cv::Mat::ones(3, 3, CV_8U);
  cv::dilate(dist, dist, kernel1);
  // cv::imshow("Dilated Distance map", dist);

  cv::waitKey(0);

  cv::Scalar ICDOpenCV = cv::mean(dist);
  cout << "ICD computed by openCV: " << ICDOpenCV[0] << endl;
  
  // Compute average intercapillary distance
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
  
  cout << "ICD: " << ICD/count << " with " << count << " pixels counted out of: " << dims[0]*dims[1] << "." << endl;
  
  writer->SetFileName("DistanceMap.vti");
  writer->SetInputData(distanceMap);
  writer->Write();
  
  
  return EXIT_SUCCESS;
}
