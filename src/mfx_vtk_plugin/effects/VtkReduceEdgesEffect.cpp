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
#include "VtkReduceEdgesEffect.h"

const char *VtkReduceEdgesEffect::GetName() {
    return "Reduce edges";
}

OfxStatus
VtkReduceEdgesEffect::vtkDescribe(OfxParamSetHandle parameters, MfxInputDef &input_mesh, MfxInputDef &output_mesh) {
    AddParam(PARAM_MAXIMUM_LENGTH, 1.0).Range(0, 1e6).Label("Maximum length");
    return kOfxStatOK;
}

OfxStatus VtkReduceEdgesEffect::vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) {

    auto maximum_length = GetParam<double>(PARAM_MAXIMUM_LENGTH).GetValue();

    return vtkCook_inner(input_polydata, output_polydata, maximum_length);
}

OfxStatus VtkReduceEdgesEffect::vtkCook_inner(vtkPolyData *input_polydata, vtkPolyData *output_polydata,
                                              double maximum_length) {
    auto vtk_input_lines = input_polydata->GetLines();
    auto vtk_output_lines = vtkSmartPointer<vtkCellArray>::New();
    auto temp_polydata = vtkSmartPointer<vtkPolyData>::New();
    temp_polydata->SetPoints(input_polydata->GetPoints());
    temp_polydata->SetLines(vtk_output_lines);

    vtk_input_lines->ConvertTo32BitStorage();
    int removed_edge_count = 0;

    for (int i = 0; i < vtk_input_lines->GetNumberOfCells(); i++) {
        int j = vtk_input_lines->GetOffsetsArray32()->GetValue(i);
        int n = vtk_input_lines->GetOffsetsArray32()->GetValue(i+1) - j;

        if (n != 2) {
            continue; // not VTK_LINE
        }

        int p1 = vtk_input_lines->GetConnectivityArray32()->GetValue(j);
        int p2 = vtk_input_lines->GetConnectivityArray32()->GetValue(j+1);

        double a[3], b[3];
        input_polydata->GetPoints()->GetPoint(p1, a);
        input_polydata->GetPoints()->GetPoint(p2, b);

        double distance_sq = (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]);

        if (distance_sq < maximum_length*maximum_length) {
            vtk_output_lines->InsertNextCell({p1, p2});
        } else {
            removed_edge_count++;
        }
    }

    if (removed_edge_count) {
        // get rid of unused points
        auto clean_polydata_filter = vtkSmartPointer<vtkCleanPolyData>::New();
        clean_polydata_filter->SetInputData(temp_polydata);
        clean_polydata_filter->SetConvertStripsToPolys(false);
        clean_polydata_filter->SetConvertPolysToLines(false);
        clean_polydata_filter->SetConvertLinesToPoints(false);
        clean_polydata_filter->SetPointMerging(false);
        clean_polydata_filter->Update();
        output_polydata->ShallowCopy(clean_polydata_filter->GetOutput());
    } else {
        output_polydata->ShallowCopy(temp_polydata);
    }

    return kOfxStatOK;
}
