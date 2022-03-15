#include<vtkSmartPointer.h>
#include<vtkImageData.h>
#include<vtkPolyData.h>
#include<vtkXMLPolyDataReader.h>
#include<vtkBox.h>
#include<vtkClipPolyData.h>
#include<vtkXMLImageDataWriter.h>

// Read the source file.
vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
reader->SetFileName("/mnt/c/Users/rhernand/Desktop/MyCodes/Virtual-Angiogenesis/perfused_geometries/Baseline.vtp");
reader->Update();  // Needed because of GetScalarRange
vtkPolyData *polyData = reader->GetOutput();

// Clipping
vtkSmartPointer<vtkBox> boundingBox = vtkSmartPointer<vtkBox>::New();
boundingBox->SetBounds(-0.15, 0.15,
		       -0.15, 0.15,
		       -0.15, 0.15);

vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
clipper->SetInputData(polyData);
clipper->SetClipFunction(boundingBox);
clipper->SetInsideOut(True);
clipper->Update();
polyData = clipper->GetOutput();

vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();
bounds = polyData->GetBounds();
int *dims {304, 304, 1};
whiteImage->SetDimensions(dims)
whiteImage->SetScalarTypeToUnsignedChar();
whiteImage->SetNumberOfScalarComponents(1);
whiteImage->AllocateScalars();

dims = whiteImage->GetDimensions();
int z = 0;
unsigned char *pixels = static_cast<unsigned char *>(whiteImage->GetScalarPointer());
for(unsigned int x = 0; x < dims[0]; x++)
  {
    for(unsigned int y = 0; y < dims[1]; y++)
      {
	pxImageData[y * dims[0] + x] = 1;   
      }
  }

vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
writer->SetFileName('Test.vti');
writer->SetInputData(whiteImage);
writer->Write();



