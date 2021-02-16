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

OfxStatus VtkEffect::Describe(OfxMeshEffectHandle descriptor) {
    printf("== VtkEffect::Describe (%s @ %p)\n", GetName(), this);
    input_definitions.clear();

    auto input_mesh = vtkAddInput(kOfxMeshMainInput);
    auto output_mesh = vtkAddInput(kOfxMeshMainOutput, true);

    OfxParamSetHandle parameters;
    meshEffectSuite->getParamSet(descriptor, &parameters);

    auto status = vtkDescribe(parameters, *input_mesh, *output_mesh);
    if (status != kOfxStatOK) {
        printf("==/ VtkEffect::Describe (failed)\n");
        return status;
    }

    // run definitions through the actual OpenMeshEffect C++ API
    for (auto &input_def : input_definitions) {
        auto native_input_def = AddInput(input_def->name);

        if (input_def->label) {
            native_input_def.Label(input_def->label);
        }

        if (!input_def->is_output) {
            native_input_def.RequestGeometry(input_def->request_geometry);
            native_input_def.RequestTransform(input_def->request_transform);
        }

        for (auto &attr_def : input_def->requested_attributes) {
            native_input_def.RequestAttribute(attr_def.attachment,
                                              attr_def.name,
                                              attr_def.componentCount,
                                              attr_def.type,
                                              attr_def.semantic,
                                              attr_def.mandatory);
        }
    }

    printf("==/ VtkEffect::Describe\n");
    return status;
}

OfxStatus VtkEffect::Cook(OfxMeshEffectHandle instance) {
    printf("== VtkEffect::Cook (%s @ %p)\n", GetName(), this);

    auto dt = [](std::chrono::system_clock::time_point t1, std::chrono::system_clock::time_point t2) -> int {
        return std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    };

    std::vector<VtkEffectInput> vtk_inputs;
    std::vector<MfxMesh> mfx_input_meshes_only; // only input meshes, no OfxMainOutputMesh
    vtk_inputs.reserve(input_definitions.size());

    VtkEffectInput *vtk_main_input = nullptr, *vtk_main_output = nullptr;

    int t_mfx_prologue = 0;
    int t_mfx_epilogue = 0;
    int t_vtk_prologue = 0;
    int t_vtk_epilogue = 0;

    for (auto &input_def : input_definitions) {
        vtk_inputs.push_back(VtkEffectInput{
            .data = nullptr,
            .definition = input_def.get()
        });
    }

    auto t_cook_start = std::chrono::system_clock::now();

    // prepare input, MFX -> VTK
    for (auto &vtk_input : vtk_inputs) {
        if (vtk_input.definition->name == kOfxMeshMainInput) {
            vtk_main_input = &vtk_input;
        } else if (vtk_input.definition->name == kOfxMeshMainOutput) {
            vtk_main_output = &vtk_input;
        }

        if (vtk_input.definition->is_output) {
            // prepare empty VTK output
            vtk_input.data = vtkSmartPointer<vtkPolyData>::New();
        } else {
            // convert input from MFX
            auto t_mfx_start = std::chrono::system_clock::now();
            MfxMesh input_mesh = GetInput(vtk_input.definition->name).GetMesh();
            auto t_mfx_end = std::chrono::system_clock::now();
            mfx_input_meshes_only.push_back(input_mesh);

            if (input_mesh.IsValid()) {
                double *transform_tmp = nullptr;
                input_mesh.FetchTransform(&transform_tmp);

                if (transform_tmp != nullptr) {
                    for (int i = 0; i < 16; i++) {
                        vtk_input.transform[i] = transform_tmp[i];
                    }
                }

                auto t_vtk_start = std::chrono::system_clock::now();
                mfx_mesh_to_vtkpolydata(vtk_input, input_mesh);
                auto t_vtk_end = std::chrono::system_clock::now();

                t_vtk_prologue += dt(t_vtk_start, t_vtk_end);
            } else {
                printf("VtkEffect::Cook - input '%s' is not valid, skipping MFX -> VTK conversion\n", vtk_input.definition->name);
            }

            t_mfx_prologue += dt(t_mfx_start, t_mfx_end);
        }
    }

    // do the job, VTK
    auto t_cook_before_vtk_cook = std::chrono::system_clock::now();
    OfxStatus cook_status = vtkCook(*vtk_main_input, *vtk_main_output, vtk_inputs);
    auto t_cook_after_vtk_cook = std::chrono::system_clock::now();

    if (cook_status != kOfxStatOK) {
        for (auto &input_mesh : mfx_input_meshes_only) {
            if (input_mesh.IsValid()) {
                input_mesh.Release();
            }
        }
      printf("==/ VtkEffect::Cook (failed)\n");
        return cook_status;
    }

    // prepare output and release it, VTK -> MFX
    // TODO support multiple output meshes, if it's ever relevant
    {
        auto t_mfx_start = std::chrono::system_clock::now();
        MfxMesh output_mesh = GetInput(kOfxMeshMainOutput).GetMesh();
        auto t_mfx_end = std::chrono::system_clock::now();
        t_mfx_epilogue += dt(t_mfx_start, t_mfx_end);

        auto t_vtk_start = std::chrono::system_clock::now();
        vtkpolydata_to_mfx_mesh(*vtk_main_output, output_mesh);
        auto t_vtk_end = std::chrono::system_clock::now();
        t_vtk_epilogue += dt(t_vtk_start, t_vtk_end);

        // "release" the MFX mesh - this will convert output data to host
        // note: avoid calling GetMesh() for output twice.
        t_mfx_start = std::chrono::system_clock::now();
        output_mesh.Release();
        t_mfx_end = std::chrono::system_clock::now();
        t_mfx_epilogue += dt(t_mfx_start, t_mfx_end);
    }

    // release inputs, MFX
    // (be sure to release output first, to prevent issues with attribute forwarding in MFX)
    auto t_mfx_start = std::chrono::system_clock::now();
    for (auto &input_mesh : mfx_input_meshes_only) {
        // make sure not to call GetMesh() multiple times for given input mesh;
        // it triggers converting on MFX Blender host side again...
        // this is why we cache input MfxMeshes
        if (input_mesh.IsValid()) {
            input_mesh.Release();
        }
    }

    auto t_mfx_end = std::chrono::system_clock::now();
    t_mfx_epilogue += dt(t_mfx_start, t_mfx_end);

    auto t_cook_after_mfx_epilogue = std::chrono::system_clock::now();

    // end
    int t_brutto = dt(t_cook_start, t_cook_after_mfx_epilogue);
    int t_netto = dt(t_cook_before_vtk_cook, t_cook_after_vtk_cook);

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

VtkEffectInputDef* VtkEffect::vtkAddInput(const char *name, bool is_output) {
    VtkEffectInputDef *ptr = new VtkEffectInputDef(name, is_output);
    input_definitions.emplace_back(ptr);
    return ptr;
}

VtkEffectInput *VtkEffect::vtkFindInput(std::vector<VtkEffectInput> &extra_inputs, const char *name) {
    VtkEffectInput *input_ptr = nullptr;

    for (int i = 0; i < extra_inputs.size(); i++) {
        if (strcmp(extra_inputs[i].definition->name, name) == 0 && extra_inputs[i].data) {
            input_ptr = &extra_inputs[i];
        }
    }

    return input_ptr;
}
