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

#include <vtkAppendPolyData.h>
#include <vtkPolyDataPointSampler.h>
#include <vtkTriangleFilter.h>

#include "VtkSamplePointsSurfaceEffect.h"

const char *VtkSamplePointsSurfaceEffect::GetName() {
    return "Sample points (surface)";
}

OfxStatus VtkSamplePointsSurfaceEffect::vtkDescribe(OfxParamSetHandle parameters, VtkEffectInputDef &input_mesh,
                                                    VtkEffectInputDef &output_mesh) {
    AddParam(PARAM_DISTANCE, 0.1).Range(1e-6, 1e6).Label("Distance");
    // AddParam(PARAM_DISTRIBUTE_UNIFORMLY, true).Label("Distribute points uniformly");
    // AddParam(PARAM_GENERATE_VERTEX_POINTS, true).Label("Generate vertex points");
    AddParam(PARAM_GENERATE_EDGE_POINTS, true).Label("Sample edges");
    AddParam(PARAM_GENERATE_INTERIOR_POINTS, true).Label("Sample faces");
    // AddParam(PARAM_INTERPOLATE_POINT_DATA, false).Label("Interpolate point data");
    return kOfxStatOK;
}

OfxStatus VtkSamplePointsSurfaceEffect::vtkCook(VtkEffectInput &main_input, VtkEffectInput &main_output, std::vector<VtkEffectInput> &extra_inputs) {
    auto distance = GetParam<double>(PARAM_DISTANCE).GetValue();
    // auto distribute_uniformly = GetParam<bool>(PARAM_DISTRIBUTE_UNIFORMLY).GetValue();
    auto generate_vertex_points = true; // GetParam<bool>(PARAM_GENERATE_VERTEX_POINTS).GetValue(); // TODO false crashes VTK 9.0.1, why?
    auto generate_edge_points = GetParam<bool>(PARAM_GENERATE_EDGE_POINTS).GetValue();
    auto generate_interior_points = GetParam<bool>(PARAM_GENERATE_INTERIOR_POINTS).GetValue();
    // bool interpolate_point_data = GetParam<bool>(PARAM_INTERPOLATE_POINT_DATA).GetValue();

    return vtkCook_inner(main_input.data, main_output.data, distance, generate_vertex_points,
                         generate_edge_points, generate_interior_points);
}

OfxStatus
VtkSamplePointsSurfaceEffect::vtkCook_inner(vtkPolyData *input_polydata, vtkPolyData *output_polydata, double distance,
                                            bool generate_vertex_points, bool generate_edge_points,
                                            bool generate_interior_points) {
    auto append_poly_data = vtkSmartPointer<vtkAppendPolyData>::New();

    auto vertex_edge_sampler = vtkSmartPointer<vtkPolyDataPointSampler>::New();
    vertex_edge_sampler->SetInputData(input_polydata);
    vertex_edge_sampler->SetDistance(distance);
    vertex_edge_sampler->SetGenerateVertexPoints(generate_vertex_points);
    vertex_edge_sampler->SetGenerateEdgePoints(generate_edge_points);
    vertex_edge_sampler->SetGenerateInteriorPoints(false);
    // vertex_edge_sampler->SetInterpolatePointData(interpolate_point_data);
    vertex_edge_sampler->Update();

    append_poly_data->AddInputData(vertex_edge_sampler->GetOutput());

    if (generate_interior_points) {
        // to handle non-convex polygons correctly, we need to triangulate first; fixes #2
        auto triangle_filter = vtkSmartPointer<vtkTriangleFilter>::New();
        triangle_filter->SetInputData(input_polydata);

        auto face_sampler = vtkSmartPointer<vtkPolyDataPointSampler>::New();
        face_sampler->SetInputConnection(triangle_filter->GetOutputPort());
        face_sampler->SetDistance(distance);
        face_sampler->SetGenerateVertexPoints(false);
        face_sampler->SetGenerateEdgePoints(false);
        face_sampler->SetGenerateInteriorPoints(true);
        // face_sampler->SetInterpolatePointData(interpolate_point_data);

        face_sampler->Update();
        append_poly_data->AddInputData(face_sampler->GetOutput());
    }

    append_poly_data->Update();

    auto filter_output = append_poly_data->GetOutput();
    output_polydata->ShallowCopy(filter_output);
    return kOfxStatOK;
}
