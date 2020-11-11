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

#include "VtkEffect.h"
#include "VtkEffectUtils.h"

#include <chrono>
#include <vtkXMLPolyDataWriter.h>

OfxStatus VtkEffect::Describe(OfxMeshEffectHandle descriptor) {
    auto input_mesh = AddInput(kOfxMeshMainInput);
    auto output_mesh = AddInput(kOfxMeshMainOutput);

    OfxParamSetHandle parameters;
    meshEffectSuite->getParamSet(descriptor, &parameters);

    return vtkDescribe(parameters, input_mesh, output_mesh);
}

OfxStatus VtkEffect::Cook(OfxMeshEffectHandle instance) {
    printf("== VtkEffect::Cook\n");

    // prepare input, OFX side
    auto t_cook_start = std::chrono::system_clock::now();
    MfxMesh input_mesh = GetInput(kOfxMeshMainInput).GetMesh();
    auto t_cook_after_mfx_prologue = std::chrono::system_clock::now();

    /*
    printf("input transform matrix is:\n");
    MfxMeshProps mesh_props;
    input_mesh.FetchProperties(mesh_props);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            printf("\t%+.3f", mesh_props.transformMatrix[i][j]);
        }
        printf("\n");
    }
    */

    // prepare input, VTK side
    auto vtk_input_polydata = mfx_mesh_to_vtkpolydata(input_mesh);  // TODO handle attributes
    auto vtk_output_polydata = vtkSmartPointer<vtkPolyData>::New();
    auto t_cook_after_vtk_prologue = std::chrono::system_clock::now();

    // debug, save input mesh
//    {
//        // vtk_input_polydata->Print(std::cout);
//        auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
//        writer->SetInputData(vtk_input_polydata);
//        writer->SetFileName("/tmp/test-input.vtp");
//        writer->Write();
//    }

    // do the job, VTK side
    auto t_cook_before_vtk_cook = std::chrono::system_clock::now();
    OfxStatus cook_status = vtkCook(vtk_input_polydata, vtk_output_polydata);
    auto t_cook_after_vtk_cook = std::chrono::system_clock::now();

    if (cook_status != kOfxStatOK) {
        input_mesh.Release();
        return cook_status;
    }

    // debug, save output mesh
//    {
//        // vtk_output_polydata->Print(std::cout);
//        auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
//        writer->SetInputData(vtk_output_polydata);
//        writer->SetFileName("/tmp/test-output.vtp");
//        writer->Write();
//    }

    // prepare output, VTK side
    auto t_cook_before_vtk_epilogue = std::chrono::system_clock::now();
    MfxMesh output_mesh = GetInput(kOfxMeshMainOutput).GetMesh();
    vtkpolydata_to_mfx_mesh(output_mesh, vtk_output_polydata);  // TODO handle attributes
    auto t_cook_after_vtk_epilogue = std::chrono::system_clock::now();

    // prepare output, OFX side
    output_mesh.Release();
    input_mesh.Release();
    auto t_cook_after_mfx_epilogue = std::chrono::system_clock::now();

    // end
    auto dt = [](std::chrono::system_clock::time_point t1, std::chrono::system_clock::time_point t2) -> int {
        return std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    };

    int t_brutto = dt(t_cook_start, t_cook_after_mfx_epilogue);
    int t_netto = dt(t_cook_before_vtk_cook, t_cook_after_vtk_cook);
    int t_mfx_prologue = dt(t_cook_start, t_cook_after_mfx_prologue);
    int t_mfx_epilogue = dt(t_cook_after_vtk_epilogue, t_cook_after_mfx_epilogue);
    int t_vtk_prologue = dt(t_cook_after_mfx_prologue, t_cook_after_vtk_prologue);
    int t_vtk_epilogue = dt(t_cook_before_vtk_epilogue, t_cook_after_vtk_epilogue);

    printf("\n\tVtkEffect cooked in %d ms (+ %d ms = %d/%d ms VTK, %d/%d ms MFX prologue/epilogue)\n\n",
           t_netto, t_brutto-t_netto, t_vtk_prologue, t_vtk_epilogue, t_mfx_prologue, t_mfx_epilogue);
    printf("==/ VtkEffect::Cook\n");
    return kOfxStatOK;
}

OfxStatus VtkEffect::IsIdentity(OfxMeshEffectHandle descriptor) {
    OfxParamSetHandle parameters;
    meshEffectSuite->getParamSet(descriptor, &parameters);

    if (vtkIsIdentity(parameters)) {
        printf("VtkEffect::IsIdentity - yes\n");
        return kOfxStatOK;
    } else {
        printf("VtkEffect::IsIdentity - no\n");
        return kOfxStatReplyDefault;
    }
}

bool VtkEffect::vtkIsIdentity(OfxParamSetHandle parameters) {
    return false;
}
