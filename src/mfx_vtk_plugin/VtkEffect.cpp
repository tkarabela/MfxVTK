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
#include <vtkCellArrayIterator.h>

static void strided_copy(void *dest_ptr, const void *src_ptr, int size, int count, int dest_stride, int src_stride) {
    assert(dest_ptr != nullptr);
    assert(src_ptr != nullptr);
    assert(size > 0);
    assert(count >= 0);
    assert(dest_stride >= size);
    assert(src_stride >= size);

    unsigned char *dest = static_cast<unsigned char *>(dest_ptr);
    const unsigned char *src = static_cast<const unsigned char *>(src_ptr);

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
    printf("== VtkEffect::Cook\n");

    auto t_cook_start = std::chrono::system_clock::now();
    MfxMesh input = GetInput(kOfxMeshMainInput).GetMesh();

    MfxMeshProps inputProps;
    input.FetchProperties(inputProps);

    MfxAttributeProps pointPos, vertPoint, faceLen;
    input.GetPointAttribute(kOfxMeshAttribPointPosition).FetchProperties(pointPos);
    input.GetVertexAttribute(kOfxMeshAttribVertexPoint).FetchProperties(vertPoint);
    input.GetFaceAttribute(kOfxMeshAttribFaceCounts).FetchProperties(faceLen);

    printf("MFX input has %d points, %d vertices, %d faces\n",
           inputProps.pointCount, inputProps.vertexCount, inputProps.faceCount);

    auto t_cook_after_mfx_prologue = std::chrono::system_clock::now();

    // ------------------------------------------------------------------------
    // Create vtkPolyData from MFX mesh
    // ------------------------------------------------------------------------
    auto vtk_input_polydata = vtkSmartPointer<vtkPolyData>::New();
    auto vtk_input_points = vtkSmartPointer<vtkPoints>::New();
    auto vtk_input_polys = vtkSmartPointer<vtkCellArray>::New();
    auto vtk_input_lines = vtkSmartPointer<vtkCellArray>::New();
    // handle isolated vertices - do we want to make vtkVertex:es out of them? (perhaps not)

    vtk_input_polydata->SetPoints(vtk_input_points);
    vtk_input_polydata->SetLines(vtk_input_lines);
    vtk_input_polydata->SetPolys(vtk_input_polys);

    // make VTK use the same datatypes as MFX
    vtk_input_points->SetDataTypeToFloat();
    vtk_input_lines->Use32BitStorage();
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

    // pre-allocation hints
    if (inputProps.vertexCount >= 3*inputProps.faceCount) {
        // looks like polygonal mesh, expect 100% polys
        vtk_input_polys->AllocateExact(inputProps.faceCount, inputProps.vertexCount);
    } else {
        // at least some faces have < 3 vertices, expect 100% lines
        vtk_input_lines->AllocateExact(inputProps.faceCount, inputProps.vertexCount);
    }

    // TODO take advantage of kOfxMeshPropNoLooseEdge, kOfxMeshPropConstantFaceCount
    int loose_edge_count = 0;

    int vertex_idx = 0;
    for (int face_idx = 0; face_idx < inputProps.faceCount; face_idx++) {
        int vertcount = (-1 == inputProps.constantFaceCount) ?
                        (*(int*)(faceLen.data + face_idx*faceLen.stride)) :
                        inputProps.constantFaceCount;

        vtkCellArray *cell_arr = nullptr;

        if (vertcount == 2) {
            loose_edge_count++;
            cell_arr = vtk_input_lines;
        } else if (vertcount >= 3) {
            cell_arr = vtk_input_polys;
        }

        if (cell_arr) {
            cell_arr->InsertNextCell(vertcount);
        }

        for (int i = 0; i < vertcount; i++, vertex_idx++) {
            int vertex_arr = *((int*)(&vertPoint.data[vertex_idx * vertPoint.stride]));
            if (cell_arr) {
                cell_arr->InsertCellPoint(vertex_arr);
            }
        }
    }

    printf("MFX input has %d loose edges\n", loose_edge_count);

    // TODO handle attributes

    vtk_input_polydata->Print(std::cout);

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

    auto t_cook_after_vtk_prologue = std::chrono::system_clock::now();
    OfxStatus cook_status = vtkCook(vtk_input_polydata, vtk_output_polydata);
    auto t_cook_after_vtk_cook = std::chrono::system_clock::now();

    if (cook_status != kOfxStatOK) {
        input.Release();
        return cook_status;
    }

//    {
//        auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
//        writer->SetInputData(vtk_output_polydata);
//        writer->SetFileName("mfx_vtk_output.vtp");
//        writer->Write();
//        printf("Dumped output to mfx_vtk_output.vtp\n");
//    }

    // ------------------------------------------------------------------------
    // Create MFX mesh from vtkPolyData
    //
    // Mapping is as follows:
    // vertex    <- VTK: vtkVertex, vtkPolyVertex
    //           -> MFX: point (without vertex or face)
    // edge      <- VTK: vtkLine, vtkPolyLine
    //           -> MFX: 2-vertex face(s)
    // face      <- VTK: vtkTriangle, vtkQuad, vtkPolygon
    //           -> MFX: n-vertex face
    // ------------------------------------------------------------------------

    auto vtk_output_verts = vtk_output_polydata->GetVerts();
    auto vtk_output_lines = vtk_output_polydata->GetLines();
    auto vtk_output_polys = vtk_output_polydata->GetPolys();
    auto vtk_output_points = vtk_output_polydata->GetPoints();

    int output_point_count = 0;
    int output_vertex_count = 0;
    int output_face_count = 0;
    int output_loose_edge_count = 0;

    if (vtk_output_points != nullptr) {
        output_point_count += vtk_output_points->GetNumberOfPoints();
    }

    if (vtk_output_verts != nullptr) {
        // no contribution to vertex or face count
    }
    if (vtk_output_lines != nullptr) {
        output_loose_edge_count = (
                vtk_output_lines->GetConnectivityArray()->GetNumberOfValues() - 1 // edge count if it was all one line
                - vtk_output_lines->GetNumberOfCells() + 1 // subtract edges between distinct lines
                );

        output_vertex_count += 2*output_loose_edge_count;
        output_face_count += output_loose_edge_count;
    }
    if (vtk_output_polys != nullptr) {
        output_vertex_count += vtk_output_polys->GetNumberOfConnectivityIds();
        output_face_count += vtk_output_polys->GetNumberOfCells();
    }

    // Allocate output MFX mesh
    int output_no_loose_edge = (output_loose_edge_count > 0) ? 0 : 1;
    int output_constant_face_count = -1;  // TODO take advantage of this
    MfxMesh output = GetInput(kOfxMeshMainOutput).GetMesh();
    output.Allocate(output_point_count, output_vertex_count, output_face_count,
                    output_no_loose_edge, output_constant_face_count); // not a hint, this must be exact

    // Get output mesh data
    output.GetPointAttribute(kOfxMeshAttribPointPosition).FetchProperties(pointPos);
    output.GetVertexAttribute(kOfxMeshAttribVertexPoint).FetchProperties(vertPoint);
    output.GetFaceAttribute(kOfxMeshAttribFaceCounts).FetchProperties(faceLen);

    printf("MFX output has %d points, %d vertices, %d faces\n",
           output_point_count, output_vertex_count, output_face_count);

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

    // write points
    if (vtk_output_points != nullptr) {
        for (int point_idx = 0; point_idx < output_point_count; point_idx++) {
            float *p = ((float *) (&pointPos.data[point_idx * pointPos.stride]));
            double p_in[3];
            vtk_output_points->GetPoint(point_idx, p_in);
            p[0] = (float)p_in[0];
            p[1] = (float)p_in[1];
            p[2] = (float)p_in[2];
        }
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

    vertex_idx = 0;
    int face_idx = 0;

    // write line vertices+faces
    if (vtk_output_lines != nullptr && vtk_output_lines->GetNumberOfCells() > 0) {
        printf("Writing lines to MFX output mesh\n");
        // slow path, accessing buffers directly would be faster
        auto iter = vtk::TakeSmartPointer(vtk_output_lines->NewIterator());
        for (iter->GoToFirstCell(); !iter->IsDoneWithTraversal(); iter->GoToNextCell()) {
            vtkIdType cellSize;
            const vtkIdType* cellPoints;
            iter->GetCurrentCell(cellSize, cellPoints);

            int* mfx_output_count = (int*)(&faceLen.data[face_idx * faceLen.stride]);
            int* mfx_output_vertex = (int*)(&vertPoint.data[vertex_idx * vertPoint.stride]);

            for (int i = 0; i < cellSize - 1; i++) {
                *mfx_output_count = 2;
                face_idx++;

                mfx_output_vertex[0] = cellPoints[i];
                mfx_output_vertex[1] = cellPoints[i+1];
                vertex_idx += 2;
            }
        }
    }

    // write poly vertices+faces
    if (vtk_output_polys != nullptr && vtk_output_polys->GetNumberOfCells() > 0) {
        printf("Writing polys to MFX output mesh\n");
        vtk_output_polys->ConvertTo32BitStorage();
        int poly_face_idx = 0;
        for (int i = 0; i < output_face_count; i++, face_idx++, poly_face_idx++) {
            int* mfx_output_count = ((int*)(&faceLen.data[face_idx * faceLen.stride]));

            int vtk_offset_start = vtk_output_polys->GetOffsetsArray32()->GetValue(poly_face_idx);
            int vtk_offset_end = vtk_output_polys->GetOffsetsArray32()->GetValue(poly_face_idx+1);
            int vtk_output_count = vtk_offset_end - vtk_offset_start;

            *mfx_output_count = vtk_output_count;

            // begin cell
            // XXX check if its polygon or something else
            for (int j = 0; j < vtk_output_count; j++, vertex_idx++) {
                int* mfx_output_vertex = ((int*)(&vertPoint.data[vertex_idx * vertPoint.stride]));

                int vtk_output_vertex = vtk_output_polys->GetConnectivityArray32()->GetValue(vtk_offset_start + j);
                *mfx_output_vertex = vtk_output_vertex;
            }
            // end cell
        }
    }

    // TODO handle attributes

    auto t_cook_after_vtk_epilogue = std::chrono::system_clock::now();

    // Release meshes
    input.Release();
    output.Release();

    auto t_cook_after_mfx_epilogue = std::chrono::system_clock::now();

    auto dt = [](std::chrono::system_clock::time_point t1, std::chrono::system_clock::time_point t2) -> int {
        return std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    };

    int t_brutto = dt(t_cook_start, t_cook_after_mfx_epilogue);
    int t_netto = dt(t_cook_after_vtk_prologue, t_cook_after_vtk_cook);
    int t_mfx_prologue = dt(t_cook_start, t_cook_after_mfx_prologue);
    int t_mfx_epilogue = dt(t_cook_after_vtk_epilogue, t_cook_after_mfx_epilogue);
    int t_vtk_prologue = dt(t_cook_after_mfx_prologue, t_cook_after_vtk_prologue);
    int t_vtk_epilogue = dt(t_cook_after_vtk_cook, t_cook_after_vtk_epilogue);

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
