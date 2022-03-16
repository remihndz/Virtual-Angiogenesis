#include<IntercapillaryDistance.h>

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
  double spacing[3] = {(bounds[1]-bounds[0])/n, (bounds[3]-bounds[2])/n, 1};


  vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
  clipper->SetInputData(tree);
  clipper->SetClipFunction(boundingBox);
  clipper->SetInsideOut(true);
  clipper->Update();
  
  vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
  int *dims = dimensions;
  dims[0] = n;
  dims[1] = n;
  dims[2] = 1;
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
  
  vtkSmartPointer<vtkPoints> points = clipper->GetOutput()->GetPoints();
  vtkSmartPointer<vtkCellArray> lines = clipper->GetOutput()->GetLines();
  vtkIdType *indices;
  vtkIdType numberOfPoints;
  unsigned int lineCount = 0;
  for (lines->InitTraversal(); lines->GetNextCell(numberOfPoints, indices); lineCount++)
    {
      // Set pixels on the line as foreground
      double* x0, x1;
      points->GetPoint(indices[0], x0);
      points->GetPoint(indices[1], x1);
      double dt = 1e-5;
      for (double t = 0
      
  
// To compute intercapillary distance from a .vtp file
double IntercapillaryDistance(const char *);
