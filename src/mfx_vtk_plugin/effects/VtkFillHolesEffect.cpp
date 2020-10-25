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

#include <vtkFillHolesFilter.h>

#include "mfx_vtk_utils.h"
#include "VtkFillHolesEffect.h"

const char *VtkFillHolesEffect::GetName() {
    return "Fill holes";
}

OfxStatus VtkFillHolesEffect::vtkDescribe(OfxParamSetHandle parameters) {
    AddParam(PARAM_HOLE_SIZE, 1.0).Range(0, 1e6).Label("Maximum hole size");
    return kOfxStatOK;
}

bool VtkFillHolesEffect::vtkIsIdentity(OfxParamSetHandle parameters) {
    double hole_size = GetParam<double>(PARAM_HOLE_SIZE).GetValue();
    return !is_positive_double(hole_size);
}

OfxStatus VtkFillHolesEffect::vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) {
    auto hole_size = GetParam<double>(PARAM_HOLE_SIZE).GetValue();
    return vtkCook_inner(input_polydata, output_polydata, hole_size);
}

OfxStatus
VtkFillHolesEffect::vtkCook_inner(vtkPolyData *input_polydata, vtkPolyData *output_polydata, double hole_size) {
    auto filter = vtkSmartPointer<vtkFillHolesFilter>::New();
    filter->SetInputData(input_polydata);

    filter->SetHoleSize(hole_size);

    filter->Update();

    auto filter_output = filter->GetOutput();
    output_polydata->ShallowCopy(filter_output);
    return kOfxStatOK;
}
