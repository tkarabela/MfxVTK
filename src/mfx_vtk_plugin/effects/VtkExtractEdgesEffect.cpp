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

#include <vtkFeatureEdges.h>
#include <vtkTriangleFilter.h>
#include <vtkCleanPolyData.h>
#include "VtkExtractEdgesEffect.h"

const char *VtkExtractEdgesEffect::GetName() {
    return "Feature edges";
}

OfxStatus VtkExtractEdgesEffect::vtkDescribe(OfxParamSetHandle parameters) {
    AddParam(PARAM_FEATURE_ANGLE, 30.0).Range(0, 180.0).Label("Feature angle");
    AddParam(PARAM_FEATURE_EDGES, true).Label("Extract feature edges");
    AddParam(PARAM_BOUNDARY_EDGES, false).Label("Extract boundary edges");
    AddParam(PARAM_NONMANIFOLD_EDGES, false).Label("Extract non-manifold edges");
    AddParam(PARAM_MANIFOLD_EDGES, false).Label("Extract manifold edges");
    return kOfxStatOK;
}

OfxStatus VtkExtractEdgesEffect::vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) {
    auto feature_angle = GetParam<double>(PARAM_FEATURE_ANGLE).GetValue();
    auto extract_feature_edges = GetParam<bool>(PARAM_FEATURE_EDGES).GetValue();
    auto extract_boundary_edges = GetParam<bool>(PARAM_BOUNDARY_EDGES).GetValue();
    auto extract_nonmanifold_edges = GetParam<bool>(PARAM_NONMANIFOLD_EDGES).GetValue();
    auto extract_manifold_edges = GetParam<bool>(PARAM_MANIFOLD_EDGES).GetValue();

    return vtkCook_inner(input_polydata, output_polydata, feature_angle, extract_feature_edges,
                         extract_boundary_edges, extract_nonmanifold_edges, extract_manifold_edges);
}

OfxStatus
VtkExtractEdgesEffect::vtkCook_inner(vtkPolyData *input_polydata, vtkPolyData *output_polydata, double feature_angle,
                                     bool extract_feature_edges, bool extract_boundary_edges,
                                     bool extract_nonmanifold_edges, bool extract_manifold_edges) {
    auto feature_edges_filter = vtkSmartPointer<vtkFeatureEdges>::New();

    if (extract_feature_edges && !extract_manifold_edges) {
        // Extract_feature_edges only works with triangles in vtkFeatureEdges, so we have to triangulate.
        // Since all feature edges are manifold edges, we can skip triangulation if extract_manifold_edges is on.
        // Triangulation does not affect non-manifold and boundary edges.
        auto triangle_filter = vtkSmartPointer<vtkTriangleFilter>::New();
        triangle_filter->SetInputData(input_polydata);
        triangle_filter->Update();
        feature_edges_filter->SetInputConnection(triangle_filter->GetOutputPort());
    } else {
        feature_edges_filter->SetInputData(input_polydata);
    }

    printf("feature_angle = %g\n", feature_angle);
    feature_edges_filter->SetFeatureAngle(feature_angle);
    feature_edges_filter->SetFeatureEdges(extract_feature_edges);
    feature_edges_filter->SetBoundaryEdges(extract_boundary_edges);
    feature_edges_filter->SetNonManifoldEdges(extract_nonmanifold_edges);
    feature_edges_filter->SetManifoldEdges(extract_manifold_edges);
    feature_edges_filter->SetColoring(false);

    // get rid of unused points
    auto clean_polydata_filter = vtkSmartPointer<vtkCleanPolyData>::New();
    clean_polydata_filter->SetInputConnection(feature_edges_filter->GetOutputPort());
    clean_polydata_filter->SetConvertStripsToPolys(false);
    clean_polydata_filter->SetConvertPolysToLines(false);
    clean_polydata_filter->SetConvertLinesToPoints(false);
    clean_polydata_filter->SetPointMerging(false);

    clean_polydata_filter->Update();

    auto filter_output = clean_polydata_filter->GetOutput();
    output_polydata->ShallowCopy(filter_output);
    return kOfxStatOK;
}
