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

#include "VtkEffectUtils.h"

#include <vtkXMLPolyDataWriter.h>
#include <vtkCellArrayIterator.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <cassert>
#include <chrono>

template <typename T, int num_components>
static void strided_copy_serial(void *dest_ptr, const void *src_ptr, int count, int dest_stride, int src_stride) {
    for (int i = 0; i < count; i++) {
        const T* src = reinterpret_cast<const T*>(reinterpret_cast<const char*>(src_ptr) + i*src_stride);
        T* dest = reinterpret_cast<T*>(reinterpret_cast<char*>(dest_ptr) + i*dest_stride);
        for (int j = 0; j < num_components; j++) {
            dest[j] = src[j];
        }
    }
}

template <typename T, int num_components>
static void strided_copy_parallel(void *dest_ptr, const void *src_ptr, int count, int dest_stride, int src_stride) {
    #pragma omp parallel for schedule(static, 1000) default(none) shared(count, dest_stride, src_stride, dest_ptr, src_ptr)
    for (int i = 0; i < count; i++) {
        const T* src = reinterpret_cast<const T*>(reinterpret_cast<const char*>(src_ptr) + i*src_stride);
        T* dest = reinterpret_cast<T*>(reinterpret_cast<char*>(dest_ptr) + i*dest_stride);
        for (int j = 0; j < num_components; j++) {
            dest[j] = src[j];
        }
    }
}

template <typename T, int num_components>
static void strided_copy(void *dest_ptr, const void *src_ptr, int count, int dest_stride, int src_stride) {
    constexpr int size = num_components*sizeof(T);

    static_assert(size > 0, "size of one element must be positive");
    assert(dest_ptr != nullptr);
    assert(src_ptr != nullptr);
    assert(count >= 0);
    assert(dest_stride >= size);
    assert(src_stride >= size);

    if (count > 5000) {
        strided_copy_parallel<T, num_components>(dest_ptr, src_ptr, count, dest_stride, src_stride);
    } else {
        strided_copy_serial<T, num_components>(dest_ptr, src_ptr, count, dest_stride, src_stride);
    }
}

vtkSmartPointer<vtkPolyData> mfx_mesh_to_vtkpolydata_generic(MfxMesh &input_mesh) {
    MfxMeshProps inputProps;
    input_mesh.FetchProperties(inputProps);

    MfxAttributeProps pointPos, vertPoint, faceLen;
    input_mesh.GetPointAttribute(kOfxMeshAttribPointPosition).FetchProperties(pointPos);
    input_mesh.GetVertexAttribute(kOfxMeshAttribVertexPoint).FetchProperties(vertPoint);
    input_mesh.GetFaceAttribute(kOfxMeshAttribFaceCounts).FetchProperties(faceLen);

    printf("MFX input has %d points, %d vertices, %d faces\n",
           inputProps.pointCount, inputProps.vertexCount, inputProps.faceCount);

    // ------------------------------------------------------------------------
    // Create vtkPolyData from MFX mesh
    // ------------------------------------------------------------------------
    auto vtk_input_polydata = vtkSmartPointer<vtkPolyData>::New();
    auto vtk_input_points = vtkSmartPointer<vtkPoints>::New();
    auto vtk_input_polys = vtkSmartPointer<vtkCellArray>::New();
    auto vtk_input_lines = vtkSmartPointer<vtkCellArray>::New();

    vtk_input_polydata->SetPoints(vtk_input_points);
    vtk_input_polydata->SetLines(vtk_input_lines);
    vtk_input_polydata->SetPolys(vtk_input_polys);

    // make VTK use the same datatypes as MFX
    vtk_input_points->SetDataTypeToFloat();
    vtk_input_lines->Use32BitStorage();
    vtk_input_polys->Use32BitStorage();

    // copy points
    vtk_input_points->SetNumberOfPoints(inputProps.pointCount);
    strided_copy<float, 3>(vtk_input_points->GetVoidPointer(0),
                           pointPos.data,
                           inputProps.pointCount,
                           sizeof(float[3]),
                           pointPos.stride);

    // copy lines + polys
    int loose_edge_count = 0;
    int vertex_idx = 0;
    for (int face_idx = 0; face_idx < inputProps.faceCount; face_idx++) {
        int vertcount = (-1 == inputProps.constantFaceCount) ?
                        (*(int *) (faceLen.data + face_idx * faceLen.stride)) :
                        inputProps.constantFaceCount;

        vtkCellArray *cell_arr = nullptr;

        if (vertcount == 2) {
            loose_edge_count++;
            cell_arr = vtk_input_lines;
            // printf("face[%d] is edge\n", face_idx);
        } else if (vertcount >= 3) {
            cell_arr = vtk_input_polys;
            // printf("face[%d] is poly %d\n", face_idx, vertcount);
        } else {
            // printf("face[%d] is invalid! %d\n", face_idx, vertcount);
        }

        if (cell_arr) {
            cell_arr->InsertNextCell(vertcount);
        }

        for (int i = 0; i < vertcount; i++, vertex_idx++) {
            int vertex_arr = *((int *) (&vertPoint.data[vertex_idx * vertPoint.stride]));
            if (cell_arr) {
                cell_arr->InsertCellPoint(vertex_arr);
            }
        }
    }
    printf("MFX input has %d loose edges\n", loose_edge_count);

    return vtk_input_polydata;
}

// pre: no loose edges
vtkSmartPointer<vtkPolyData> mfx_mesh_to_vtkpolydata_polys(MfxMesh &input_mesh) {
    MfxMeshProps inputProps;
    input_mesh.FetchProperties(inputProps);

    MfxAttributeProps pointPos, vertPoint, faceLen;
    input_mesh.GetPointAttribute(kOfxMeshAttribPointPosition).FetchProperties(pointPos);
    input_mesh.GetVertexAttribute(kOfxMeshAttribVertexPoint).FetchProperties(vertPoint);
    input_mesh.GetFaceAttribute(kOfxMeshAttribFaceCounts).FetchProperties(faceLen);

    printf("MFX input has %d points, %d vertices, %d faces\n",
           inputProps.pointCount, inputProps.vertexCount, inputProps.faceCount);

    // ------------------------------------------------------------------------
    // Create vtkPolyData from MFX mesh
    // ------------------------------------------------------------------------
    auto vtk_input_polydata = vtkSmartPointer<vtkPolyData>::New();
    auto vtk_input_points = vtkSmartPointer<vtkPoints>::New();
    auto vtk_input_polys = vtkSmartPointer<vtkCellArray>::New();

    vtk_input_polydata->SetPoints(vtk_input_points);
    vtk_input_polydata->SetPolys(vtk_input_polys);

    // make VTK use the same datatypes as MFX
    vtk_input_points->SetDataTypeToFloat();
    vtk_input_polys->Use32BitStorage();

    // copy points
    if (inputProps.pointCount > 0) {
        vtk_input_points->SetNumberOfPoints(inputProps.pointCount);
        strided_copy<float, 3>(vtk_input_points->GetVoidPointer(0),
                               pointPos.data,
                               inputProps.pointCount,
                               sizeof(float[3]),
                               pointPos.stride);
    }

    // copy vertices
    if (inputProps.vertexCount > 0) {
        vtk_input_polys->GetConnectivityArray32()->SetNumberOfValues(inputProps.vertexCount);
        strided_copy<int, 1>(vtk_input_polys->GetConnectivityArray32()->GetVoidPointer(0),
                             vertPoint.data,
                             inputProps.vertexCount,
                             sizeof(int),
                             vertPoint.stride);
    }

    // copy faces
    if (inputProps.faceCount > 0) {
        vtk_input_polys->GetOffsetsArray32()->SetNumberOfValues(inputProps.faceCount + 1);
        int *offset_array = reinterpret_cast<int*>(vtk_input_polys->GetOffsetsArray32()->GetVoidPointer(0));
        int vertex_sum = 0;
        if (inputProps.constantFaceCount == -1) {
            // varying face counts
            for (int i = 0; i < inputProps.faceCount; i++) {
                offset_array[i] = vertex_sum;
                vertex_sum += *reinterpret_cast<int*>(faceLen.data + i*faceLen.stride);
            }
        } else {
            for (int i = 0; i < inputProps.faceCount; i++) {
                offset_array[i] = vertex_sum;
                vertex_sum += inputProps.constantFaceCount;
            }
        }
        offset_array[inputProps.faceCount] = vertex_sum; // sentinel value
    }

    // copy UV and color attributes
    // TODO switch this to more general approach when available; implement for other mesh types as well
    // note we cannot forward here since we have per-vertex data from OFX but VTK only has per-point data
    for (int k = 0; k < 4; k++) {
        char name[32];
        sprintf(name, "color%d", k);

        if (input_mesh.HasVertexAttribute(name)) {
            MfxAttributeProps attr;
            input_mesh.GetVertexAttribute(name).FetchProperties(attr);
            assert(attr.type == MfxAttributeType::UByte);

            auto array = vtkSmartPointer<vtkUnsignedCharArray>::New();
            array->SetNumberOfComponents(attr.componentCount);
            array->SetNumberOfTuples(vtk_input_polydata->GetNumberOfPoints());
            array->FillValue(0);
            array->SetName(name);
            vtk_input_polydata->GetPointData()->AddArray(array);

            for (int i = 0; i < inputProps.vertexCount; i++) {
                int p = *(int*)(vertPoint.data + i*vertPoint.stride);
                char *values = attr.data + i*attr.stride;
                for (int j = 0; j < attr.componentCount; j++) {
                    array->SetTypedComponent(p, j, values[j]);
                }
            }
            printf("MfxVTK - read array %s\n", name);
        }
    }

    for (int k = 0; k < 4; k++) {
        char name[32];
        sprintf(name, "uv%d", k);

        if (input_mesh.HasVertexAttribute(name)) {
            MfxAttributeProps attr;
            input_mesh.GetVertexAttribute(name).FetchProperties(attr);
            assert(attr.type == MfxAttributeType::Float);

            auto array = vtkSmartPointer<vtkFloatArray>::New();
            array->SetNumberOfComponents(attr.componentCount);
            array->SetNumberOfTuples(vtk_input_polydata->GetNumberOfPoints());
            array->FillValue(0);
            array->SetName(name);
            vtk_input_polydata->GetPointData()->AddArray(array);

            for (int i = 0; i < inputProps.vertexCount; i++) {
                int p = *(int*)(vertPoint.data + i*vertPoint.stride);
                float *values = (float*)(attr.data + i*attr.stride);
                for (int j = 0; j < attr.componentCount; j++) {
                    array->SetTypedComponent(p, j, values[j]);
                }
            }
            printf("MfxVTK - read array %s\n", name);
        }
    }

    return vtk_input_polydata;
}

vtkSmartPointer<vtkPolyData> mfx_mesh_to_vtkpolydata(MfxMesh &input_mesh) {
    MfxMeshProps inputProps;
    input_mesh.FetchProperties(inputProps);

    if (inputProps.noLooseEdge) {
        return mfx_mesh_to_vtkpolydata_polys(input_mesh);
    } else {
        return mfx_mesh_to_vtkpolydata_generic(input_mesh);
    }
}

// pre: no lines/polys
static void vtkpolydata_to_mfx_mesh_pointcloud(MfxMesh &output_mesh, vtkPolyData *vtk_output_polydata) {
    auto attrib_point_position = output_mesh.GetPointAttribute(kOfxMeshAttribPointPosition);
    auto attrib_vertex_point = output_mesh.GetVertexAttribute(kOfxMeshAttribVertexPoint);
    auto attrib_face_counts = output_mesh.GetFaceAttribute(kOfxMeshAttribFaceCounts);

    MfxAttributeProps attrib_point_position_props, attrib_vertex_point_props, attrib_face_counts_props;
    attrib_point_position.FetchProperties(attrib_point_position_props);
    attrib_vertex_point.FetchProperties(attrib_vertex_point_props);
    attrib_face_counts.FetchProperties(attrib_face_counts_props);

    int point_count = vtk_output_polydata->GetNumberOfPoints();
    int vertex_count = 0;
    int face_count = 0;
    int no_loose_edge = 1;
    int constant_face_count = 3; // dummy value

    vtk_output_polydata->GetPoints()->SetDataTypeToFloat();
    attrib_point_position_props.isOwner = false;
    attrib_point_position_props.data = reinterpret_cast<char*>(vtk_output_polydata->GetPoints()->GetVoidPointer(0));
    attrib_point_position_props.stride = 3*sizeof(float);

    attrib_vertex_point_props.isOwner = false;
    attrib_vertex_point_props.data = nullptr;

    attrib_face_counts_props.isOwner = false;
    attrib_face_counts_props.data = nullptr;

    attrib_point_position.SetProperties(attrib_point_position_props);
    attrib_vertex_point.SetProperties(attrib_vertex_point_props);
    attrib_face_counts.SetProperties(attrib_face_counts_props);

    output_mesh.Allocate(point_count, vertex_count, face_count, no_loose_edge, constant_face_count);
}

// pre: no polys, only vtkLines, no vtkPolyLine
static void vtkpolydata_to_mfx_mesh_wireframe(MfxMesh &output_mesh, vtkPolyData *vtk_output_polydata) {
    auto attrib_point_position = output_mesh.GetPointAttribute(kOfxMeshAttribPointPosition);
    auto attrib_vertex_point = output_mesh.GetVertexAttribute(kOfxMeshAttribVertexPoint);
    auto attrib_face_counts = output_mesh.GetFaceAttribute(kOfxMeshAttribFaceCounts);

    MfxAttributeProps attrib_point_position_props, attrib_vertex_point_props, attrib_face_counts_props;
    attrib_point_position.FetchProperties(attrib_point_position_props);
    attrib_vertex_point.FetchProperties(attrib_vertex_point_props);
    attrib_face_counts.FetchProperties(attrib_face_counts_props);

    int vtk_line_count = vtk_output_polydata->GetNumberOfLines();
    int point_count = vtk_output_polydata->GetNumberOfPoints();
    int vertex_count = 2*vtk_line_count;
    int face_count = vtk_line_count;
    int no_loose_edge = 0;
    int constant_face_count = 2;

    vtk_output_polydata->GetPoints()->SetDataTypeToFloat();
    attrib_point_position_props.isOwner = false;
    attrib_point_position_props.data = reinterpret_cast<char*>(vtk_output_polydata->GetPoints()->GetVoidPointer(0));
    attrib_point_position_props.stride = 3*sizeof(float);

    vtk_output_polydata->GetLines()->ConvertTo32BitStorage();
    attrib_vertex_point_props.isOwner = false;
    attrib_vertex_point_props.data = reinterpret_cast<char*>(vtk_output_polydata->GetLines()->GetConnectivityArray32()->GetVoidPointer(0));
    attrib_vertex_point_props.stride = sizeof(int);

    attrib_face_counts_props.isOwner = false;
    attrib_face_counts_props.data = nullptr;

    attrib_point_position.SetProperties(attrib_point_position_props);
    attrib_vertex_point.SetProperties(attrib_vertex_point_props);
    attrib_face_counts.SetProperties(attrib_face_counts_props);

    output_mesh.Allocate(point_count, vertex_count, face_count, no_loose_edge, constant_face_count);
}

// pre: only polys, no lines
static void vtkpolydata_to_mfx_mesh_poly(MfxMesh &output_mesh, vtkPolyData *vtk_output_polydata) {
    auto attrib_point_position = output_mesh.GetPointAttribute(kOfxMeshAttribPointPosition);
    auto attrib_vertex_point = output_mesh.GetVertexAttribute(kOfxMeshAttribVertexPoint);
    auto attrib_face_counts = output_mesh.GetFaceAttribute(kOfxMeshAttribFaceCounts);

    MfxAttributeProps attrib_point_position_props, attrib_vertex_point_props, attrib_face_counts_props;
    attrib_point_position.FetchProperties(attrib_point_position_props);
    attrib_vertex_point.FetchProperties(attrib_vertex_point_props);
    attrib_face_counts.FetchProperties(attrib_face_counts_props);

    auto t0 = std::chrono::system_clock::now();
    int point_count = vtk_output_polydata->GetNumberOfPoints();
    int vertex_count = vtk_output_polydata->GetPolys()->GetConnectivityArray()->GetNumberOfValues();
    int face_count = vtk_output_polydata->GetPolys()->GetNumberOfCells();
    int no_loose_edge = 1;
    int constant_face_count = (vtk_output_polydata->GetPolys()->GetMaxCellSize() == 3) ? 3 : -1;
    auto t1 = std::chrono::system_clock::now();

    printf("vtkpolydata_to_mfx_mesh forwarding points\n");
    vtk_output_polydata->GetPoints()->SetDataTypeToFloat();
    attrib_point_position_props.isOwner = false;
    attrib_point_position_props.data = reinterpret_cast<char*>(vtk_output_polydata->GetPoints()->GetVoidPointer(0));
    attrib_point_position_props.stride = 3*sizeof(float);
    auto t2 = std::chrono::system_clock::now();

    printf("vtkpolydata_to_mfx_mesh forwarding vertices\n");
    vtk_output_polydata->GetPolys()->ConvertTo32BitStorage();
    attrib_vertex_point_props.isOwner = false;
    attrib_vertex_point_props.data = reinterpret_cast<char*>(vtk_output_polydata->GetPolys()->GetConnectivityArray32()->GetVoidPointer(0));
    attrib_vertex_point_props.stride = sizeof(int);
    auto t3 = std::chrono::system_clock::now();

    if (constant_face_count == -1) {
        printf("vtkpolydata_to_mfx_mesh copying face counts\n");
        attrib_face_counts_props.isOwner = true; // request buffer, VTK offset array cannot be used unfortunately
    } else {
        printf("vtkpolydata_to_mfx_mesh forwarding face counts\n");
        attrib_face_counts_props.isOwner = false;
        attrib_face_counts_props.data = nullptr;
    }

    // handle vertex colors, UVs
    // TODO do this more generally, for more mesh types
    for (int k = 0; k < 4; k++) {
        char name[32];
        sprintf(name, "color%d", k);
        auto array = vtk_output_polydata->GetPointData()->GetArray(name);
        if (array != nullptr) {
            printf("vtkpolydata_to_mfx_mesh copying attribute %s\n", name);
            output_mesh.AddVertexAttribute(name, array->GetNumberOfComponents(), kOfxMeshAttribTypeUByte);
        }
    }
    for (int k = 0; k < 4; k++) {
        char name[32];
        sprintf(name, "uv%d", k);
        auto array = vtk_output_polydata->GetPointData()->GetArray(name);
        if (array != nullptr) {
            printf("vtkpolydata_to_mfx_mesh copying attribute %s\n", name);
            output_mesh.AddVertexAttribute(name, array->GetNumberOfComponents(), kOfxMeshAttribTypeFloat);
        }
    }

    attrib_point_position.SetProperties(attrib_point_position_props);
    attrib_vertex_point.SetProperties(attrib_vertex_point_props);
    attrib_face_counts.SetProperties(attrib_face_counts_props);

    printf("vtkpolydata_to_mfx_mesh allocate - point_count %d, vertex_count %d, face_count %d, no_loose_edge %d, constant_face_count %d\n",
           point_count, vertex_count, face_count, no_loose_edge, constant_face_count);
    output_mesh.Allocate(point_count, vertex_count, face_count, no_loose_edge, constant_face_count);
    auto t4 = std::chrono::system_clock::now();

    if (constant_face_count == -1) {
        attrib_face_counts.FetchProperties(attrib_face_counts_props);

        int *offset_array = reinterpret_cast<int*>(vtk_output_polydata->GetPolys()->GetOffsetsArray32()->GetVoidPointer(0));

        for (int i = 0; i < face_count; i++) {
            int *count = reinterpret_cast<int*>(attrib_face_counts_props.data + i*attrib_face_counts_props.stride);
            *count = offset_array[i+1] - offset_array[i];
        }
    }
    auto t5 = std::chrono::system_clock::now();

    // handle vertex colors, UVs
    // TODO do this more generally, for more mesh types
    // note we cannot forward here since we have per-point data but OFX expects per-vertex data
    for (int k = 0; k < 4; k++) {
        char name[32];
        sprintf(name, "color%d", k);
        auto array = vtk_output_polydata->GetPointData()->GetArray(name);
        if (array != nullptr) {
            MfxAttributeProps attr;
            output_mesh.GetVertexAttribute(name).FetchProperties(attr);

            for (int i = 0; i < vertex_count; i++) {
                int p = *(int*)(attrib_vertex_point_props.data + i*attrib_vertex_point_props.stride);
                for (int j = 0; j < array->GetNumberOfComponents(); j++) {
                    attr.data[i*attr.stride + j] = array->GetComponent(p, j); // TODO use GetTypedComponent here
                }
            }
            printf("MfxVTK - wrote array %s\n", name);
        }
    }
    for (int k = 0; k < 4; k++) {
        char name[32];
        sprintf(name, "uv%d", k);
        auto array = vtk_output_polydata->GetPointData()->GetArray(name);
        if (array != nullptr) {
            MfxAttributeProps attr;
            output_mesh.GetVertexAttribute(name).FetchProperties(attr);

            for (int i = 0; i < vertex_count; i++) {
                int p = *(int*)(attrib_vertex_point_props.data + i*attrib_vertex_point_props.stride);
                float *dest = (float*)(attr.data + i*attr.stride);
                for (int j = 0; j < array->GetNumberOfComponents(); j++) {
                    dest[j] = array->GetComponent(p, j);
                }
            }
            printf("MfxVTK - wrote array %s\n", name);
        }
    }

    auto dt = [](std::chrono::system_clock::time_point t1, std::chrono::system_clock::time_point t2) -> int {
        return std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    };

    printf("\n\tvtkpolydata_to_mfx_mesh_poly timings [ms] dt1=%d, dt2=%d, dt3=%d, dt4=%d, dt5=%d\n",
           dt(t0, t1), dt(t1, t2), dt(t2, t3), dt(t3, t4), dt(t4, t5));
}

// pre: nothing
static void vtkpolydata_to_mfx_mesh_generic(MfxMesh &output_mesh, vtkPolyData *vtk_output_polydata) {
    auto vtk_output_verts = vtk_output_polydata->GetVerts();
    auto vtk_output_lines = vtk_output_polydata->GetLines();
    auto vtk_output_polys = vtk_output_polydata->GetPolys();
    auto vtk_output_points = vtk_output_polydata->GetPoints();

    int output_point_count = 0;
    int output_vertex_count = 0;
    int output_face_count = 0;
    int output_loose_edge_count = 0;
    int output_poly_count = 0;

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
        output_poly_count += vtk_output_polys->GetNumberOfCells();
    }

    // Allocate output MFX mesh
    int output_no_loose_edge = (output_loose_edge_count > 0) ? 0 : 1;
    int output_constant_face_count = (output_poly_count == 0 && output_loose_edge_count > 0) ? 2 : -1;  // TODO take further advantage of this
    output_mesh.Allocate(output_point_count, output_vertex_count, output_face_count,
                         output_no_loose_edge, output_constant_face_count); // not a hint, this must be exact

    // Get output mesh data
    MfxAttributeProps pointPos, vertPoint, faceLen;
    output_mesh.GetPointAttribute(kOfxMeshAttribPointPosition).FetchProperties(pointPos);
    output_mesh.GetVertexAttribute(kOfxMeshAttribVertexPoint).FetchProperties(vertPoint);
    output_mesh.GetFaceAttribute(kOfxMeshAttribFaceCounts).FetchProperties(faceLen);

    printf("MFX output has %d points, %d vertices, %d faces\n",
           output_point_count, output_vertex_count, output_face_count);

    // copy points
    if (vtk_output_points != nullptr) {
        strided_copy<float, 3>(pointPos.data,
                               vtk_output_points->GetVoidPointer(0),
                               output_point_count,
                               pointPos.stride,
                               sizeof(float[3]));
    }

    // copy cells
    int vertex_idx = 0;
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
}

void vtkpolydata_to_mfx_mesh(MfxMesh &output_mesh, vtkPolyData *vtk_output_polydata) {
    auto vtk_output_lines = vtk_output_polydata->GetLines();
    auto vtk_output_polys = vtk_output_polydata->GetPolys();
    bool has_lines = (vtk_output_lines == nullptr || vtk_output_lines->GetNumberOfCells() > 0);
    bool has_polylines = (has_lines && vtk_output_lines->GetMaxCellSize() > 2);
    bool has_polys = (vtk_output_polys == nullptr || vtk_output_polys->GetNumberOfCells() > 0);

    if (has_lines) {
        if (!has_polys && !has_polylines) {
            vtkpolydata_to_mfx_mesh_wireframe(output_mesh, vtk_output_polydata);
        } else {
            vtkpolydata_to_mfx_mesh_generic(output_mesh, vtk_output_polydata);
        }
    } else {
        // no lines
        if (has_polys) {
            vtkpolydata_to_mfx_mesh_poly(output_mesh, vtk_output_polydata);
        } else {
            vtkpolydata_to_mfx_mesh_pointcloud(output_mesh, vtk_output_polydata);
        }
    }
}
