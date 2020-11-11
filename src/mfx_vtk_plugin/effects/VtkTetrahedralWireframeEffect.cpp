/*
MfxVTK Open Mesh Effect plug-in
Copyright (c) 2020 Tomas Karabela

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#include <vtkCleanPolyData.h>
#include <vtkDelaunay3D.h>
#include <vtkExtractEdges.h>
#include <vtkUnstructuredGrid.h>

#include "VtkReduceEdgesEffect.h"
#include "VtkTetrahedralWireframeEffect.h"

const char *VtkTetrahedralWireframeEffect::GetName() {
    return "Tetrahedral wireframe";
}

OfxStatus VtkTetrahedralWireframeEffect::vtkDescribe(OfxParamSetHandle parameters, MfxInputDef &input_mesh,
                                                     MfxInputDef &output_mesh) {
    AddParam(PARAM_MAXIMUM_EDGE_LENGTH, 1.0).Range(0, 1e6).Label("Maximum edge length");
    return kOfxStatOK;
}

OfxStatus VtkTetrahedralWireframeEffect::vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) {
    auto maximum_edge_length = GetParam<double>(PARAM_MAXIMUM_EDGE_LENGTH).GetValue();
    return vtkCook_inner(input_polydata, output_polydata, maximum_edge_length);
}

OfxStatus VtkTetrahedralWireframeEffect::vtkCook_inner(vtkPolyData *input_polydata, vtkPolyData *output_polydata,
                                                       double maximum_edge_length) {
    auto delaunay3d_filter = vtkSmartPointer<vtkDelaunay3D>::New();
    delaunay3d_filter->SetInputData(input_polydata);
    delaunay3d_filter->SetAlpha(0);
    delaunay3d_filter->Update();

    auto extract_edges = vtkSmartPointer<vtkExtractEdges>::New();
    extract_edges->SetInputData(delaunay3d_filter->GetOutput());
    extract_edges->Update();

    auto raw_edge_polydata = extract_edges->GetOutput();

    return VtkReduceEdgesEffect::vtkCook_inner(raw_edge_polydata, output_polydata, maximum_edge_length);
}
