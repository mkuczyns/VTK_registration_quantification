cmake_minimum_required(VERSION 2.8)

PROJECT(vtkRegistrationQuantification)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(vtkRegistrationQuantification MACOSX_BUNDLE vtkRegistrationQuantification.cxx)

if(VTK_LIBRARIES)
  target_link_libraries(vtkRegistrationQuantification ${VTK_LIBRARIES})
else()
  target_link_libraries(vtkRegistrationQuantification vtkHybrid vtkWidgets)
endif()

