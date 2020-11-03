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

set(SRC
  VtkEffect.cpp
  VtkEffect.h
        mfx_vtk_utils.h
  effects/VtkReduceEdgesEffect.cpp effects/VtkReduceEdgesEffect.h effects/VtkTetrahedralWireframeEffect.cpp effects/VtkTetrahedralWireframeEffect.h effects/VtkSamplePointsVolumeEffect.cpp effects/VtkSamplePointsVolumeEffect.h effects/VtkMakeTubesEffect.cpp effects/VtkMakeTubesEffect.h effects/VtkFillHolesEffect.cpp effects/VtkFillHolesEffect.h effects/VtkSamplePointsSurfaceEffect.cpp effects/VtkSamplePointsSurfaceEffect.h effects/VtkDecimateEffect.cpp effects/VtkDecimateEffect.h effects/VtkSmoothEffect.cpp effects/VtkSmoothEffect.h effects/VtkSurfaceDistanceEffect.cpp effects/VtkSurfaceDistanceEffect.h VtkEffectUtils.cpp VtkEffectUtils.h)

set(LIB
  openmesheffect_openfx
  CppPluginSupport
)

# vtk
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

add_library(mfx_vtk_plugin SHARED mfx_vtk_plugin.cpp ${SRC} effects/VtkExtractEdgesEffect.cpp effects/VtkExtractEdgesEffect.h)
target_include_directories(mfx_vtk_plugin PRIVATE ${INC})
target_link_libraries(mfx_vtk_plugin PRIVATE ${LIB} ${VTK_LIBRARIES})
set_target_properties(mfx_vtk_plugin PROPERTIES SUFFIX ".ofx")

add_library(mfx_vtk_plugin_extra SHARED mfx_vtk_plugin_extra.cpp ${SRC})
target_include_directories(mfx_vtk_plugin_extra PRIVATE ${INC})
target_link_libraries(mfx_vtk_plugin_extra PRIVATE ${LIB} ${VTK_LIBRARIES})
set_target_properties(mfx_vtk_plugin_extra PROPERTIES SUFFIX ".ofx")

vtk_module_autoinit(
        TARGETS mfx_vtk_plugin mfx_vtk_plugin_extra
        MODULES ${VTK_LIBRARIES}
)

# openmp
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(mfx_vtk_plugin PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(mfx_vtk_plugin_extra PUBLIC OpenMP::OpenMP_CXX)
endif()
