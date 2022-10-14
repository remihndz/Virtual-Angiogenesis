cmake_minimum_required(VERSION 2.6)

project(WholeRetinaSVP)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(WholeRetinaSVP WholeRetinaSVP.cxx)
if(VTK_LIBRARIES)
  target_link_libraries(WholeRetinaSVP ${VTK_LIBRARIES})
else()
  target_link_libraries(WholeRetinaSVP vtkHybrid)
  endif()
