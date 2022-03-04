import sys, glob
from sys import platform as _platform
import vtkmodules.all as vtk

if _platform=='linux':
    file_names = sys.argv[1:]
else:
    file_names = glob.glob(sys.argv[1])

for file_name in file_names:
    file_name_base = file_name[:-4]

    colors = vtk.vtkNamedColors()

    # Read the source file.
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(file_name)
    reader.Update()  # Needed because of GetScalarRange
    polyData = reader.GetOutput()


    # Clipping
    boundingBox = vtk.vtkBox()
    boundingBox.SetBounds(-0.15, 0.15,
                          -0.15, 0.15,
                          -0.15, 0.15)
    
    clipper = vtk.vtkClipPolyData()
    clipper.SetInputData(polyData)
    clipper.SetClipFunction(boundingBox)
    clipper.SetInsideOut(True)
    clipper.Update()
    
    polyData = clipper.GetOutput()
    
    # scalar_range = polyData.GetScalarRange()
    
    # Create the mapper that corresponds the objects of the vtk.vtk file
    # into graphics elements
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(polyData)
    mapper.ScalarVisibilityOff()

    # Create the Actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().EdgeVisibilityOn()
    actor.GetProperty().SetLineWidth(2.0)
    actor.GetProperty().SetColor(colors.GetColor3d("MistyRose"))
    
    backface = vtk.vtkProperty()
    backface.SetColor(colors.GetColor3d('Tomato'))
    actor.SetBackfaceProperty(backface)
    
    # Create the Renderer
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(colors.GetColor3d('Black'))
    
    # Create the RendererWindow
    renderer_window = vtk.vtkRenderWindow()
    renderer_window.SetSize(640, 480)
    renderer_window.AddRenderer(renderer)
    
    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput(renderer_window)
    w2if.Update()
    
    writer = vtk.vtkPNGWriter()
    writer.SetFileName(file_name_base + '.png')
    writer.SetInputConnection(w2if.GetOutputPort())
    writer.Write()

    # Create the RendererWindowInteractor and display the vtk_file
    # interactor = vtk.vtkRenderWindowInteractor()
    # interactor.SetRenderWindow(renderer_window)
    # interactor.Initialize()
    # interactor.Start()

    print('File ', file_name, 'saved as', file_name_base + '.png')