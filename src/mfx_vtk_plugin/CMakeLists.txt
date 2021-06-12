# MfxVTK Open Mesh Effect plug-in
# Copyright (c) 2020 Tomas Karabela
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

set(INC
  .
)

file(GLOB SRC_EFFECTS effects/*.cpp effects/*.h)

set(SRC
        VtkEffect.cpp
        VtkEffect.h
        VtkEffectUtils.cpp
        VtkEffectUtils.h
        VtkEffectInput.cpp
        VtkEffectInput.h
        VtkEffectInputDef.cpp
        VtkEffectInputDef.h
        mfx_vtk_utils.h
        ${SRC_EFFECTS})

# OpenMfx dependency
set(LIB
  OpenMfx::Core
  OpenMfx::Sdk
)

# Enable C++20 (required to use designated initializers)
set(CMAKE_CXX_STANDARD 20)
if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest")
endif(MSVC)

# VTK dependency, part 1
find_package(VTK COMPONENTS
        vtkCommonCore
        vtkCommonDataModel
        vtkFiltersCore
        vtkFiltersSources
        vtkFiltersModeling
        vtkFiltersGeometry
        vtkFiltersExtraction
        vtkFiltersFlowPaths
        vtkIOXML
        QUIET)
if (NOT VTK_FOUND)
    message("VTK not found: ${VTK_NOT_FOUND_MESSAGE}")
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")

# Main target (for effects in ./effects)
add_library(mfx_vtk_plugin SHARED mfx_vtk_plugin.cpp ${SRC})
target_include_directories(mfx_vtk_plugin PRIVATE ${INC})
target_link_libraries(mfx_vtk_plugin PRIVATE ${LIB} ${VTK_LIBRARIES})
set_target_properties(mfx_vtk_plugin PROPERTIES SUFFIX ".ofx")

# Extra target (for effects in mfx_vtk_plugin_extra.cpp)
add_library(mfx_vtk_plugin_extra SHARED mfx_vtk_plugin_extra.cpp ${SRC})
target_include_directories(mfx_vtk_plugin_extra PRIVATE ${INC})
target_link_libraries(mfx_vtk_plugin_extra PRIVATE ${LIB} ${VTK_LIBRARIES})
set_target_properties(mfx_vtk_plugin_extra PROPERTIES SUFFIX ".ofx")

# VTK dependency, part 2
vtk_module_autoinit(
        TARGETS mfx_vtk_plugin mfx_vtk_plugin_extra
        MODULES ${VTK_LIBRARIES}
)

# OpenMP dependency
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(mfx_vtk_plugin PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(mfx_vtk_plugin_extra PUBLIC OpenMP::OpenMP_CXX)
endif()
