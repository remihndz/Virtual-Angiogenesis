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
#include <vtkImageDilateErode3D.h>

// OpenCV libs
#include<opencv2/core.hpp>
#include<opencv2/imgproc.hpp>
#include<opencv2/highgui.hpp>
// #include<opencv2/ximgproc.hpp> // Needs specific library installation



// Functions to get pixels corresponding to x,y,z coordinates et vice-versa
void CoordinatesToPixel(double *point, int *dimensions, int* pixel);
void PixelToCoordinates(int *pixel, int *dimensions, double* point, double *bounds);

// Function to convert VTK's vtkImageData type to a openCV Mat
cv::Mat convertVtkImageDataToCVMat(const vtkSmartPointer<vtkImageData> &vtkImage);

// Function to compute intercapillary distance from a polydata structure double IntercapillaryDistance(const vtkSmartPointer<vtkPolyData> polyData); // For a FOV of 3mm and pixels=512
double IntercapillaryDistance(const vtkSmartPointer<vtkPolyData> tree, double FOV, int resolution);
double IntercapillaryDistance(const vtkSmartPointer<vtkPolyData> tree, double FOV);
// To compute intercapillary distance from a .vtp file
double IntercapillaryDistance(const char * fileNameVTP, double FOV, int resolution);
double IntercapillaryDistance(const char * fileNameVTP, double FOV);

