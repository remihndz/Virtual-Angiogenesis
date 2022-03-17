#include<IntercapillaryDistance.hpp>

// Functions to get pixels corresponding to x,y,z coordinates et vice-versa
void CoordinatesToPixel(double *point, int *dimensions, int* pixel)
{
  double dx {0.3/dimensions[0]}, dy {0.3/dimensions[1]};
  pixel[0] = floor(dimensions[0]/2.0) + int(point[0]/dx);
  pixel[1] = floor(dimensions[1]/2.0) + int(point[1]/dy);
  pixel[2] = 0;  
}

void PixelToCoordinates(int *pixel, int *dimensions, double* point, double *bounds)
{
  double dx {0.3/dimensions[0]}, dy {0.3/dimensions[1]};
  point[0] = bounds[0] + dx * pixel[0];
  point[1] = bounds[2] + dy * pixel[1];
}
;

// Function to convert VTK's vtkImageData type to a openCV Mat
cv::Mat convertVtkImageDataToCVMat(const vtkSmartPointer<vtkImageData> &vtkImage)
{
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


// Function to compute intercapillary distance from a polydata structure
double IntercapillaryDistance(const vtkSmartPointer<vtkPolyData> tree, double FOV, int resolution)
{
  double xmin {-FOV/2.0}, xmax {FOV/2.0};
  double ymin {-FOV/2.0}, ymax {FOV/2.0};
  
  unsigned int foregroundPixel {0}, backgroundPixel {255};
  
  // Clipping to FOV
  vtkSmartPointer<vtkBox> boundingBox = vtkSmartPointer<vtkBox>::New();
  boundingBox->SetBounds(xmin, xmax,
			 ymin, ymax,
			 -0.15, 0.15);
  double *bounds = boundingBox->GetBounds();
  double spacing[3] = {(bounds[1]-bounds[0])/resolution, (bounds[3]-bounds[2])/resolution, 1};

  vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
  clipper->SetInputData(tree);
  clipper->SetClipFunction(boundingBox);
  clipper->SetInsideOut(true);
  clipper->Update();
  vtkSmartPointer<vtkPolyData> clippedTree = clipper->GetOutput();
  
  vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
  int *dims = new int[3]{resolution, resolution, 1};
  // dims[0] = resolution;
  // dims[1] = resolution;
  // dims[2] = 1;
  imageData->SetDimensions(dims);
  imageData->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
  
  // Create a black canvas and print the tree in white on it
  unsigned char *pixels = static_cast<unsigned char *>(imageData->GetScalarPointer());
  
  // int extent[6];
  // imageData->GetExtent(extent);
  // cout << "Image extent: ";
  // for (const auto& e : extent)
  //   cout << e << ' ';
  // cout << endl;
  // cout << "Image dims: ";
  // for (int i = 0; i < 3; i++)
  //   cout << dims[i] << ' ';
  // cout << endl;

  for (int i = 0; i < dims[0]; i++)
    {
      for (int j = 0; j < dims[1]; j++)
	{
	  int pixel[3] = {i,j,0};
	  pixels[i * dims[0] + j] = backgroundPixel;

	  // // Alternatively, to count FAZ as 'vessels'
	  // if (pow(point[0],2) + pow(point[1],2) < pow(0.04 ,2))
	  //   pixels[i * dims[0] + j] = foregroundPixel;
	  // else
	  //  pixels[i * dims[0] + j] = backgroundPixel;   
	}
    }
 
  vtkSmartPointer<vtkPoints> points = clippedTree->GetPoints();
  vtkSmartPointer<vtkCellArray> lines = clippedTree->GetLines();
  vtkIdType *indices;
  vtkIdType numberOfPoints;
  unsigned int lineCount = 0;
  for (lines->InitTraversal(); lines->GetNextCell(numberOfPoints, indices); lineCount++)
    {
      // Set pixels on the line as foreground
      double *x0 = new double[3];
      double *x1 = new double[3];
      points->GetPoint(indices[0], x0);
      points->GetPoint(indices[1], x1);
      double dt = 1e-3;
      for (double t = 0; t<=1; t+=dt)
	{
	  double *x  = new double[3] {(1-t)*x0[0] + t*x1[0], (1-t)*x0[1] + t*x1[1], 0};
	  int *pixel = new int[3];
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
  
  // vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
  // writer->SetFileName("Test.vti");
  // writer->SetInputConnection(dilateErode->GetOutputPort());
  // writer->Write();


  // Compute the distance map using openCV
  cv::Mat bwImage = convertVtkImageDataToCVMat(dilateErode->GetOutput());
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
  // cv::waitKey(0); // Makes the plots if any

  // Compute mean ICD
  cv::Scalar ICD = cv::mean(dist);
  cout << "Intercapillary distance for the current tree: " << ICD[0] << endl;
  
  return ICD[0];
}

double IntercapillaryDistance(const vtkSmartPointer<vtkPolyData> tree, double FOV)
{
  int resolution = 1024;
  double ICD = IntercapillaryDistance(tree, FOV, resolution);
  return ICD;
}
  
// To compute intercapillary distance from a .vtp file
double IntercapillaryDistance(const char *fileNameVTP, double FOV, int resolution)
{
  vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(fileNameVTP);
  reader->Update();
  vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();
  cout << "Polydata loaded." << endl;
  double ICD;
  ICD = IntercapillaryDistance(polyData, FOV, resolution);
  return ICD;
}
 
double IntercapillaryDistance(const char *fileNameVTP, double FOV)
{
  int resolution = 1024;
  vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(fileNameVTP);
  reader->Update();
  vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();
  double ICD;
  ICD = IntercapillaryDistance(polyData, FOV, resolution);
  return ICD;
}


