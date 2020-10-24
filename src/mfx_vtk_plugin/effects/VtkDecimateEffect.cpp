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

#include <vtkTriangleFilter.h>
#include <vtkQuadricDecimation.h>

#include "VtkDecimateEffect.h"

const char *VtkDecimateEffect::GetName() {
    return "Decimate";
}

OfxStatus VtkDecimateEffect::vtkDescribe(OfxParamSetHandle parameters) {
    // to make this more akin to Blender Decimate, we use target "keep" ratio instead of "remove" ratio
    // which is used by the underlying VTK filter (SetTargetReduction)
    AddParam(PARAM_TARGET_RATIO, 1.0).Range(0.0, 1.0).Label("Target ratio");
    AddParam(PARAM_VOLUME_PRESERVATION, false).Label("Preserve volume");
    return kOfxStatOK;
}

bool VtkDecimateEffect::vtkIsIdentity(OfxParamSetHandle parameters) {
    auto target_ratio = GetParam<double>(PARAM_TARGET_RATIO).GetValue();
    return target_ratio >= 1.0;
}

OfxStatus VtkDecimateEffect::vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) {
    auto target_ratio = GetParam<double>(PARAM_TARGET_RATIO).GetValue();
    auto volume_preservation = GetParam<bool>(PARAM_VOLUME_PRESERVATION).GetValue();
    return vtkCook_inner(input_polydata, output_polydata, 1.0 - target_ratio, volume_preservation);
}

OfxStatus
VtkDecimateEffect::vtkCook_inner(vtkPolyData *input_polydata, vtkPolyData *output_polydata, double target_reduction,
                                 bool volume_preservation) {
    // vtkTriangleFilter to ensure triangle mesh on input
    auto triangle_filter = vtkSmartPointer<vtkTriangleFilter>::New();
    triangle_filter->SetInputData(input_polydata);

    // vtkQuadricDecimation for main processing
    auto decimate_filter = vtkSmartPointer<vtkQuadricDecimation>::New();
    decimate_filter->SetInputConnection(triangle_filter->GetOutputPort());
    decimate_filter->SetTargetReduction(target_reduction);
    decimate_filter->SetVolumePreservation(volume_preservation);
    // TODO the filter supports optimizing for attribute error, too, we could expose this

    decimate_filter->Update();

    auto filter_output = decimate_filter->GetOutput();
    output_polydata->ShallowCopy(filter_output);
    return kOfxStatOK;
}
