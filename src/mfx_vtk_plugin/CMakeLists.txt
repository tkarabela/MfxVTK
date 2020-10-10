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
  mfx_vtk_plugin.cpp
  VtkEffect.cpp
  VtkEffect.h)

set(LIB
  openmesheffect_openfx
  CppPluginSupport
)

add_library(mfx_vtk_plugin SHARED ${SRC})

target_include_directories(mfx_vtk_plugin PRIVATE ${INC})
target_link_libraries(mfx_vtk_plugin PRIVATE ${LIB})
set_target_properties(mfx_vtk_plugin PROPERTIES SUFFIX ".ofx")

# vtk
find_package(VTK COMPONENTS
        vtkCommonCore
        vtkCommonDataModel
        vtkFiltersCore
        vtkFiltersSources
        vtkIOXML
        QUIET)
if (NOT VTK_FOUND)
    message("VTK not found: ${VTK_NOT_FOUND_MESSAGE}")
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
target_link_libraries(mfx_vtk_plugin PRIVATE ${VTK_LIBRARIES})
vtk_module_autoinit(
        TARGETS mfx_vtk_plugin
        MODULES ${VTK_LIBRARIES}
)
