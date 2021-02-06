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

#include <vtkSmoothPolyDataFilter.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <mfx_vtk_utils.h>
#include "VtkSmoothEffect.h"

const char *VtkSmoothEffect::GetName() {
    return "Smooth";
}

OfxStatus
VtkSmoothEffect::vtkDescribe(OfxParamSetHandle parameters, VtkEffectInputDef &input_mesh, VtkEffectInputDef &output_mesh) {
    AddParam(PARAM_MODE, MODE_WINDOWED_SINC).Range(1, 2).Label("Mode"); // TODO make this enum!
    AddParam(PARAM_ITERATIONS, 20).Range(1, 1000).Label("Iterations");
    AddParam(PARAM_FACTOR, 0.1).Range(0.0, 1000.0).Label("Factor");
    AddParam(PARAM_BOUNDARY_SMOOTHING, true).Label("Boundary smoothing");
    AddParam(PARAM_FEATURE_EDGE_SMOOTHING, false).Label("Feature edge smoothing");
    AddParam(PARAM_FEATURE_ANGLE, 45.0).Range(0.001, 180.0).Label("Feature angle");
    AddParam(PARAM_EDGE_ANGLE, 15.0).Range(0.001, 180.0).Label("Edge angle");

    // TODO declare this is a deformer
    return kOfxStatOK;
}

OfxStatus VtkSmoothEffect::vtkCook(VtkEffectInput &main_input, VtkEffectInput &main_output, std::vector<VtkEffectInput> &extra_inputs) {
    auto mode = GetParam<int>(PARAM_MODE).GetValue();
    auto iterations = GetParam<int>(PARAM_ITERATIONS).GetValue();
    auto factor = GetParam<double>(PARAM_FACTOR).GetValue();
    auto boundary_smoothing = GetParam<bool>(PARAM_BOUNDARY_SMOOTHING).GetValue();
    auto feature_edge_smoothing = GetParam<bool>(PARAM_FEATURE_EDGE_SMOOTHING).GetValue();
    auto feature_angle = GetParam<double>(PARAM_FEATURE_ANGLE).GetValue();
    auto edge_angle = GetParam<double>(PARAM_EDGE_ANGLE).GetValue();

    // XXX until we have enums...
    mode = clamp(mode, 1, 2);

    if (mode == MODE_LAPLACIAN) {
        auto relaxation_factor = factor;
        return vtkCook_inner_laplacian(main_input.data, main_output.data, iterations, relaxation_factor,
                                       boundary_smoothing, feature_edge_smoothing, feature_angle, edge_angle);
    } else if (mode == MODE_WINDOWED_SINC) {
        auto passband = clamp(factor, 0.0, 2.0);
        return vtkCook_inner_windowed_sinc(main_input.data, main_output.data, iterations, passband,
                                           boundary_smoothing, feature_edge_smoothing, feature_angle, edge_angle);
    } else {
        // this should not happen
        printf("VtkSmoothEffect::vtkCook - Bad mode!\n");
        return kOfxStatErrValue;
    }
}

OfxStatus
VtkSmoothEffect::vtkCook_inner_laplacian(vtkPolyData *input_polydata, vtkPolyData *output_polydata, int iterations,
                                         double relaxation_factor, bool boundary_smoothing, bool feature_edge_smoothing,
                                         double feature_angle, double edge_angle) {
    auto filter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    filter->SetInputData(input_polydata);

    filter->SetNumberOfIterations(iterations);
    filter->SetConvergence(0.001);  // 0.1% of bounding box diagonal
    filter->SetBoundarySmoothing(boundary_smoothing);
    filter->SetRelaxationFactor(relaxation_factor);
    filter->SetFeatureEdgeSmoothing(feature_edge_smoothing);
    filter->SetFeatureAngle(feature_angle);
    filter->SetEdgeAngle(edge_angle);

    filter->Update();

    auto filter_output = filter->GetOutput();
    output_polydata->ShallowCopy(filter_output);
    return kOfxStatOK;
}

OfxStatus
VtkSmoothEffect::vtkCook_inner_windowed_sinc(vtkPolyData *input_polydata, vtkPolyData *output_polydata, int iterations,
                                             double passband, bool boundary_smoothing, bool feature_edge_smoothing,
                                             double feature_angle, double edge_angle) {
    auto filter = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
    filter->SetInputData(input_polydata);

    filter->SetNumberOfIterations(iterations);
    filter->SetPassBand(passband);
    filter->SetBoundarySmoothing(boundary_smoothing);
    filter->SetNonManifoldSmoothing(false);
    filter->SetFeatureEdgeSmoothing(feature_edge_smoothing);
    filter->SetFeatureAngle(feature_angle);
    filter->SetEdgeAngle(edge_angle);
    filter->SetNormalizeCoordinates(true);

    filter->Update();

    auto filter_output = filter->GetOutput();
    output_polydata->ShallowCopy(filter_output);
    return kOfxStatOK;
}
