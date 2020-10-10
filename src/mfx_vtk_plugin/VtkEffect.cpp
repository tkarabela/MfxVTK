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

#include <cstring>
#include <cassert>
#include <chrono>
#include <vtkXMLPolyDataWriter.h>

static void strided_copy(void *dest, void *src, int size, int count, int dest_stride, int src_stride) {
    assert(dest != nullptr);
    assert(src != nullptr);
    assert(size > 0);
    assert(count >= 0);
    assert(dest_stride >= size);
    assert(src_stride >= size);

    if (size == dest_stride && size == src_stride) {
        memcpy(dest, src, size*count);
    } else {
        for (int i = 0; i < count; i++) {
            memcpy(dest, src, size);
            dest += dest_stride;
            src += src_stride;
        }
    }
}

OfxStatus VtkEffect::Describe(OfxMeshEffectHandle descriptor) {
    AddInput(kOfxMeshMainInput);
    AddInput(kOfxMeshMainOutput);

    OfxParamSetHandle parameters;
    meshEffectSuite->getParamSet(descriptor, &parameters);

    return vtkDescribe(parameters);
}

OfxStatus VtkEffect::Cook(OfxMeshEffectHandle instance) {
    auto t_cook_start = std::chrono::system_clock::now();
    MfxMesh input = GetInput(kOfxMeshMainInput).GetMesh();

    MfxMeshProps inputProps;
    input.FetchProperties(inputProps);

    MfxAttributeProps pointPos, vertPoint, faceLen;
    input.GetPointAttribute(kOfxMeshAttribPointPosition).FetchProperties(pointPos);
    input.GetVertexAttribute(kOfxMeshAttribVertexPoint).FetchProperties(vertPoint);
    input.GetFaceAttribute(kOfxMeshAttribFaceCounts).FetchProperties(faceLen);

    // ------------------------------------------------------------------------
    // Create vtkPolyData from MFX mesh
    // TODO optimize
    // ------------------------------------------------------------------------
    auto vtk_input_polydata = vtkSmartPointer<vtkPolyData>::New();
    auto vtk_input_points = vtkSmartPointer<vtkPoints>::New();
    auto vtk_input_polys = vtkSmartPointer<vtkCellArray>::New();
    // ...and isolated vertices - not supported by MFX
    // ...and face-less edges - not supported by MFX

    vtk_input_polydata->SetPoints(vtk_input_points);
    vtk_input_polydata->SetPolys(vtk_input_polys);

    // make VTK use the same datatypes as MFX
    vtk_input_points->SetDataTypeToFloat();
    vtk_input_polys->Use32BitStorage();

    // copy points
    vtk_input_points->SetNumberOfPoints(inputProps.pointCount);
    float* vtk_input_points_buffer = static_cast<float*>(vtk_input_points->GetVoidPointer(0));

    strided_copy(vtk_input_points_buffer,
                 pointPos.data,
                 3*sizeof(float),
                 inputProps.pointCount,
                 3*sizeof(float),
                 pointPos.stride);

    // copy cells
    vtk_input_polys->AllocateExact(inputProps.faceCount, inputProps.vertexCount); // just a hint

    int vertex_idx = 0;
    for (int face_idx = 0; face_idx < inputProps.faceCount; face_idx++) {
        int vertcount = *((int*)(&faceLen.data[face_idx * faceLen.stride]));
        bool skip = false;

        if (vertcount < 3) {
            skip = true; // XXX how to handle isolated vertices/edges?
        }

        if (!skip) {
            vtk_input_polys->InsertNextCell(vertcount);
        }

        for (int i = 0; i < vertcount; i++, vertex_idx++) {
            int vertex_arr = *((int*)(&vertPoint.data[vertex_idx * vertPoint.stride]));
            if (!skip){
                vtk_input_polys->InsertCellPoint(vertex_arr);
            }
        }
    }

    // TODO handle attributes

    // ------------------------------------------------------------------------
    // Cook the vtkPolyData
    // ------------------------------------------------------------------------

//    {
//        auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
//        writer->SetInputData(vtk_input_polydata);
//        writer->SetFileName("test-1input.vtp");
//        writer->Write();
//    }

    auto vtk_output_polydata = vtkSmartPointer<vtkPolyData>::New();

    auto t_cook_inner_start = std::chrono::system_clock::now();
    OfxStatus cook_status = vtkCook(vtk_input_polydata, vtk_output_polydata);
    auto t_cook_inner_end = std::chrono::system_clock::now();

    if (cook_status != kOfxStatOK) {
        input.Release();
        return cook_status;
    }

//    {
//        auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
//        writer->SetInputData(vtk_input_polydata);
//        writer->SetFileName("test-2output.vtp");
//        writer->Write();
//    }

    // ------------------------------------------------------------------------
    // Create MFX mesh from vtkPolyData
    // TODO optimize
    // ------------------------------------------------------------------------

    auto vtk_output_polys = vtk_output_polydata->GetPolys();
    auto vtk_output_points = vtk_output_polydata->GetPoints();

    int output_point_count = vtk_output_points->GetNumberOfPoints();
    int output_vertex_count = vtk_output_polys->GetNumberOfConnectivityIds();
    int output_face_count = vtk_output_polys->GetNumberOfCells();

    // Allocate output MFX mesh
    MfxMesh output = GetInput(kOfxMeshMainOutput).GetMesh();
    output.Allocate(output_point_count, output_vertex_count, output_face_count);

    // Get output mesh data
    output.GetPointAttribute(kOfxMeshAttribPointPosition).FetchProperties(pointPos);
    output.GetVertexAttribute(kOfxMeshAttribVertexPoint).FetchProperties(vertPoint);
    output.GetFaceAttribute(kOfxMeshAttribFaceCounts).FetchProperties(faceLen);

    // copy points
    /*
    float* vtk_output_points_buffer = static_cast<float*>(vtk_output_points->GetVoidPointer(0));

    strided_copy(vtk_output_points_buffer,
                 pointPos.data,
                 3*sizeof(float),
                 output_point_count,
                 pointPos.stride,
                 3*sizeof(float));
                 */

    for (int point_idx = 0; point_idx < output_point_count; point_idx++) {
        float *p = ((float *) (&pointPos.data[point_idx * pointPos.stride]));
        double p_in[3];
        vtk_output_points->GetPoint(point_idx, p_in);
        p[0] = (float)p_in[0];
        p[1] = (float)p_in[1];
        p[2] = (float)p_in[2];
    }

    // copy cells
    /*
    vtk_output_polys->ConvertTo32BitStorage();

    int* vtk_output_connectivity_buffer = static_cast<int*>(vtk_output_polys->GetConnectivityArray()->GetVoidPointer(0));
    int* vtk_output_offsets_buffer = static_cast<int*>(vtk_output_polys->GetOffsetsArray()->GetVoidPointer(0));

    strided_copy(vtk_output_connectivity_buffer,
                 vertPoint.data,
                 sizeof(int),
                 output_vertex_count,
                 vertPoint.stride,
                 sizeof(int));

    strided_copy(vtk_output_offsets_buffer,
                 faceLen.data,
                 sizeof(int),
                 output_face_count,
                 faceLen.stride,
                 sizeof(int));
                 */

    vtk_output_polys->ConvertTo32BitStorage();
    vertex_idx = 0;
    for (int face_idx = 0; face_idx < output_face_count; face_idx++) {
        int* mfx_output_count = ((int*)(&faceLen.data[face_idx * faceLen.stride]));

        int vtk_offset_start = vtk_output_polys->GetOffsetsArray32()->GetValue(face_idx);
        int vtk_offset_end = vtk_output_polys->GetOffsetsArray32()->GetValue(face_idx+1);
        int vtk_output_count = vtk_offset_end - vtk_offset_start;

        *mfx_output_count = vtk_output_count;

        // begin cell
        // XXX check if its polygon or something else
        for (int i = 0; i < vtk_output_count; i++, vertex_idx++) {
            int* mfx_output_vertex = ((int*)(&vertPoint.data[vertex_idx * vertPoint.stride]));

            int vtk_output_vertex = vtk_output_polys->GetConnectivityArray32()->GetValue(vtk_offset_start + i);
            *mfx_output_vertex = vtk_output_vertex;
        }
        // end cell
    }

    // TODO handle attributes

    // Release meshes
    input.Release();
    output.Release();

    auto t_cook_end = std::chrono::system_clock::now();

    auto t_brutto = std::chrono::duration_cast<std::chrono::milliseconds>(t_cook_end - t_cook_start).count();
    auto t_netto = std::chrono::duration_cast<std::chrono::milliseconds>(t_cook_inner_end - t_cook_inner_start).count();

    printf("\n\tVtkEffect cooked in %d ms (+ %d ms in wrapper)\n\n", (int)t_netto, (int)(t_brutto - t_netto));

    return kOfxStatOK;
}
