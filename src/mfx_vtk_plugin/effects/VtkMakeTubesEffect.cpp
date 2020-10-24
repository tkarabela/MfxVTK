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

#include <vtkTubeFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkExtractEdges.h>

#include "VtkMakeTubesEffect.h"

const char *VtkMakeTubesEffect::GetName() {
    return "Make tubes";
}

OfxStatus VtkMakeTubesEffect::vtkDescribe(OfxParamSetHandle parameters) {
    AddParam(PARAM_RADIUS, 0.02).Range(1e-6, 1e6).Label("Radius");
    AddParam(PARAM_NUMBER_OF_SIDES, 6).Range(3, 1000).Label("Number of sides");
    AddParam(PARAM_CAPPING, true).Label("Cap ends");
    // AddParam(PARAM_TUBE_BENDER, true).Label("Use tube bender");
    // TODO texture coordinates
    return kOfxStatOK;
}

OfxStatus VtkMakeTubesEffect::vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) {
    double radius = GetParam<double>(PARAM_RADIUS).GetValue();
    int number_of_sides = GetParam<int>(PARAM_NUMBER_OF_SIDES).GetValue();
    bool capping = GetParam<bool>(PARAM_CAPPING).GetValue();
    // bool use_tube_bender = GetParam<bool>(PARAM_TUBE_BENDER).GetValue();
    return vtkCook_inner(input_polydata, output_polydata, radius, number_of_sides, capping);
}

OfxStatus VtkMakeTubesEffect::vtkCook_inner(vtkPolyData *input_polydata, vtkPolyData *output_polydata, double radius,
                                            int number_of_sides, bool capping) {
    auto append_poly_data = vtkSmartPointer<vtkAppendPolyData>::New();

    if (input_polydata->GetNumberOfPolys() > 0) {
        // vtkExtractEdges to create lines even from polygonal mesh
        auto extract_edges_filter = vtkSmartPointer<vtkExtractEdges>::New();
        extract_edges_filter->SetInputData(input_polydata);
        extract_edges_filter->Update();
        append_poly_data->AddInputData(extract_edges_filter->GetOutput());
    } else {
        append_poly_data->AddInputData(input_polydata);
    }

    // TODO incorporate optional vtkTubeBender - when it lands post VTK 9.0
    // if (use_tube_bender) {
    //    auto tube_bender_filter = vtkSmartPointer<vtkTubeBender>::New();
    //    tube_bender_filter->SetInputConnection(extract_edges_filter->GetOutputPort());
    //    tube_filter_input = tube_bender_filter->GetOutputPort();
    // }

    // vtkTubeFilter to turn lines into polygonal tubes
    auto tube_filter = vtkSmartPointer<vtkTubeFilter>::New();
    tube_filter->SetInputConnection(append_poly_data->GetOutputPort());
    tube_filter->SetRadius(radius);
    tube_filter->SetNumberOfSides(number_of_sides);
    tube_filter->SetCapping(capping);
    tube_filter->SetSidesShareVertices(true);

    // vtkTriangleFilter to convert triangle strips to polygons
    auto triangle_filter = vtkSmartPointer<vtkTriangleFilter>::New();
    triangle_filter->SetInputConnection(tube_filter->GetOutputPort());

    triangle_filter->Update();

    auto filter_output = triangle_filter->GetOutput();
    output_polydata->ShallowCopy(filter_output);
    return kOfxStatOK;
}
