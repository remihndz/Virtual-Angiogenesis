cmake_minimum_required(VERSION 3.10)

# set the project name
project(dcco_WholeRetina_SVP)
project(Macula)
project(2DRetina)


# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add all libraries needed for compilation
set(VTK_8_1_libs vtkCommonCore-8.1 vtkCommonDataModel-8.1 vtkCommonExecutionModel-8.1 vtkFiltersModeling-8.1 vtkIOCore-8.1 vtkIOLegacy-8.1 vtkIOXML-8.1 vtkIOGeometry-8.1 vtkInfovisCore-8.1 vtkFiltersGeneral-8.1 vtkFiltersCore-8.1 vtkCommonTransforms-8.1 vtkIOXMLParser-8.1)

find_path(VITA_INCLUDE_DIRS core/GeneratorData.h
          HINTS .. ENV VITA_PATH
          PATH_SUFFIXES include/vita_source )
find_path(VITA_LIBRARY_DIRS libVItA.a libVItA.so
          HINTS .. ENV VITA_PATH
          PATH_SUFFIXES lib )
find_path(VTK_INCLUDE_DIRS vtkAbstractArray.h
          HINTS .. ENV VITA_PATH /home/remi/AppsAndLibs/VItA/
          PATH_SUFFIXES vita_build/include/vtk-8.1 include/vtk-8.1)
find_path(VTK_LIBRARY_DIRS libvtkCommonCore-8.1.so
          HINTS .. ENV VITA_PATH /home/remi/AppsAndLibs/VItA/
          PATH_SUFFIXES vita_build/lib lib)

link_directories(${VITA_LIBRARY_DIRS} ${VTK_LIBRARY_DIRS})

# add the executable
add_executable(dcco_WholeRetina_SVP dcco_WholeRetina_SVP.cpp)
add_executable(Macula Macula.cpp)
add_executable(2DRetina 2DRetina.cpp)

# link libraries
find_package(OpenMP)
if (OpenMP_CXX_FOUND)
   target_link_libraries(dcco_WholeRetina_SVP PUBLIC VItA ${VTK_8_1_libs} OpenMP::OpenMP_CXX)
   target_link_libraries(Macula PUBLIC VItA ${VTK_8_1_libs} OpenMP::OpenMP_CXX)
   target_link_libraries(2DRetina PUBLIC VItA ${VTK_8_1_libs} OpenMP::OpenMP_CXX)		
else()
   target_link_libraries(dcco_WholeRetina_SVP VItA ${VTK_8_1_libs})
   target_link_libraries(Macula VItA ${VTK_8_1_libs})
   target_link_libraries(2DRetina VItA ${VTK_8_1_libs})
endif()

# add all include directories needed for compilation
target_include_directories(dcco_WholeRetina_SVP PUBLIC
                           "${VITA_INCLUDE_DIRS}"
                           "${VTK_INCLUDE_DIRS}"
                           )
target_include_directories(Macula PUBLIC
                           "${VITA_INCLUDE_DIRS}"
                           "${VTK_INCLUDE_DIRS}"
                           )
target_include_directories(2DRetina PUBLIC
                           "${VITA_INCLUDE_DIRS}"
                           "${VTK_INCLUDE_DIRS}"
                           )
