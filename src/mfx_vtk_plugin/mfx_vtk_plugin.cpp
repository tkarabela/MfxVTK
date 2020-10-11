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

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkPolyDataPointSampler.h>
#include <vtkFeatureEdges.h>
#include <vtkFillHolesFilter.h>
#include <vtkExtractEdges.h>
#include <vtkTubeFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkDecimatePro.h>

#include "ofxCore.h"
#include "ofxMeshEffect.h"

#include "PluginSupport/MfxRegister"
#include "VtkEffect.h"

// TODO implement kOfxMeshEffectActionIsIdentity for the effects

class VtkSmoothPolyDataFilterEffect : public VtkEffect {
private:
    const char *PARAM_ITERATIONS = "NumberOfIterations";
    const char *PARAM_CONVERGENCE = "Convergence";
    const char *PARAM_RELAXATION_FACTOR = "RelaxationFactor";
    const char *PARAM_BOUNDARY_SMOOTHING = "BoundarySmoothing";
    const char *PARAM_FEATURE_EDGE_SMOOTHING = "FeatureEdgeSmoothing";
    const char *PARAM_FEATURE_ANGLE = "FeatureAngle";
    const char *PARAM_EDGE_ANGLE = "EdgeAngle";

public:
    const char* GetName() override {
        return "Smooth (Laplacian)";
    }

    OfxStatus vtkDescribe(OfxParamSetHandle parameters) override {
        AddParam(PARAM_ITERATIONS, 20).Range(1, 1000).Label("Iterations");
        AddParam(PARAM_RELAXATION_FACTOR, 0.1).Range(0.0, 1000.0).Label("Relaxation factor");
        AddParam(PARAM_BOUNDARY_SMOOTHING, true).Label("Boundary smoothing");
        AddParam(PARAM_FEATURE_EDGE_SMOOTHING, false).Label("Feature edge smoothing");
        AddParam(PARAM_FEATURE_ANGLE, 45.0).Range(0.001, 180.0).Label("Feature angle");
        AddParam(PARAM_EDGE_ANGLE, 15.0).Range(0.001, 180.0).Label("Edge angle");
        AddParam(PARAM_CONVERGENCE, 0.0).Range(0.0, 1000.0).Label("Convergence");
        return kOfxStatOK;
    }

    OfxStatus vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) override {
        int iterations = GetParam<int>(PARAM_ITERATIONS).GetValue();
        double convergence = GetParam<double>(PARAM_CONVERGENCE).GetValue();
        bool boundary_smoothing = GetParam<bool>(PARAM_BOUNDARY_SMOOTHING).GetValue();
        double relaxation_factor = GetParam<double>(PARAM_RELAXATION_FACTOR).GetValue();
        bool feature_edge_smoothing = GetParam<bool>(PARAM_FEATURE_EDGE_SMOOTHING).GetValue();
        double feature_angle = GetParam<double>(PARAM_FEATURE_ANGLE).GetValue();
        double edge_angle = GetParam<double>(PARAM_EDGE_ANGLE).GetValue();

        auto filter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
        filter->SetInputData(input_polydata);

        filter->SetNumberOfIterations(iterations);
        filter->SetConvergence(convergence);
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
};

// ----------------------------------------------------------------------------

class VtkWindowedSincPolyDataFilterEffect : public VtkEffect {
private:
    const char *PARAM_ITERATIONS = "NumberOfIterations";
    const char *PARAM_PASSBAND = "PassBand";
    const char *PARAM_BOUNDARY_SMOOTHING = "BoundarySmoothing";
    const char *PARAM_NONMANIFOLD_SMOOTHING = "NonManifoldSmoothing";
    const char *PARAM_FEATURE_EDGE_SMOOTHING = "FeatureEdgeSmoothing";
    const char *PARAM_FEATURE_ANGLE = "FeatureAngle";
    const char *PARAM_EDGE_ANGLE = "EdgeAngle";
    const char *PARAM_NORMALIZE_COORDINATES = "NormalizeCoordinates";

public:
    const char* GetName() override {
        return "Smooth (windowed sinc)";
    }

    OfxStatus vtkDescribe(OfxParamSetHandle parameters) override {
        AddParam(PARAM_ITERATIONS, 20).Range(1, 1000).Label("Iterations");
        AddParam(PARAM_PASSBAND, 0.1).Range(0.001, 2.0).Label("Passband");
        AddParam(PARAM_BOUNDARY_SMOOTHING, true).Label("Boundary smoothing");
        AddParam(PARAM_NONMANIFOLD_SMOOTHING, false).Label("Non-manifold smoothing");
        AddParam(PARAM_FEATURE_EDGE_SMOOTHING, false).Label("Feature edge smoothing");
        AddParam(PARAM_FEATURE_ANGLE, 45.0).Range(0.001, 180.0).Label("Feature angle");
        AddParam(PARAM_EDGE_ANGLE, 15.0).Range(0.001, 180.0).Label("Edge angle");
        AddParam(PARAM_NORMALIZE_COORDINATES, true).Label("Normalize coordinates");
        return kOfxStatOK;
    }

    OfxStatus vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) override {
        int iterations = GetParam<int>(PARAM_ITERATIONS).GetValue();
        double passband = GetParam<double>(PARAM_PASSBAND).GetValue();
        bool boundary_smoothing = GetParam<bool>(PARAM_BOUNDARY_SMOOTHING).GetValue();
        bool nonmanifold_smoothing = GetParam<bool>(PARAM_NONMANIFOLD_SMOOTHING).GetValue();
        bool feature_edge_smoothing = GetParam<bool>(PARAM_FEATURE_EDGE_SMOOTHING).GetValue();
        double feature_angle = GetParam<double>(PARAM_FEATURE_ANGLE).GetValue();
        double edge_angle = GetParam<double>(PARAM_EDGE_ANGLE).GetValue();
        bool normalize_coordinates = GetParam<bool>(PARAM_NORMALIZE_COORDINATES).GetValue();

        auto filter = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
        filter->SetInputData(input_polydata);

        filter->SetNumberOfIterations(iterations);
        filter->SetPassBand(passband);
        filter->SetBoundarySmoothing(boundary_smoothing);
        filter->SetNonManifoldSmoothing(nonmanifold_smoothing);
        filter->SetFeatureEdgeSmoothing(feature_edge_smoothing);
        filter->SetFeatureAngle(feature_angle);
        filter->SetEdgeAngle(edge_angle);
        filter->SetNormalizeCoordinates(normalize_coordinates);

        filter->Update();

        auto filter_output = filter->GetOutput();
        output_polydata->ShallowCopy(filter_output);
        return kOfxStatOK;
    }
};

// ----------------------------------------------------------------------------

class VtkPolyDataPointSamplerEffect : public VtkEffect {
private:
    const char *PARAM_DISTANCE = "Distance";
    // const char *PARAM_USE_RANDOM_GENERATION_MODE = "UseRandomGenerationMode";
    // const char *PARAM_GENERATE_VERTEX_POINTS = "GenerateVertexPoints";
    const char *PARAM_GENERATE_EDGE_POINTS = "GenerateEdgePoints";
    const char *PARAM_GENERATE_INTERIOR_POINTS = "GenerateInteriorPoints";
    const char *PARAM_INTERPOLATE_POINT_DATA = "InterpolatePointData";

public:
    const char* GetName() override {
        return "Point sampling";
    }

    OfxStatus vtkDescribe(OfxParamSetHandle parameters) override {
        AddParam(PARAM_DISTANCE, 0.1).Range(1e-6, 1e6).Label("Distance");
        // AddParam(PARAM_USE_RANDOM_GENERATION_MODE, true).Label("Random point generation");
        // AddParam(PARAM_GENERATE_VERTEX_POINTS, true).Label("Generate vertex points");
        AddParam(PARAM_GENERATE_EDGE_POINTS, true).Label("Generate edge points");
        AddParam(PARAM_GENERATE_INTERIOR_POINTS, true).Label("Generate interior points");
        AddParam(PARAM_INTERPOLATE_POINT_DATA, false).Label("Interpolate point data");
        return kOfxStatOK;
    }

    OfxStatus vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) override {
        double distance = GetParam<double>(PARAM_DISTANCE).GetValue();
        // bool user_random_generation_mode = GetParam<bool>(PARAM_USE_RANDOM_GENERATION_MODE).GetValue();
        bool generate_vertex_points = true; // GetParam<bool>(PARAM_GENERATE_VERTEX_POINTS).GetValue(); // false crashes VTK 9.0.1, why?
        bool generate_edge_points = GetParam<bool>(PARAM_GENERATE_EDGE_POINTS).GetValue();
        bool generate_interior_points = GetParam<bool>(PARAM_GENERATE_INTERIOR_POINTS).GetValue();
        bool interpolate_point_data = GetParam<bool>(PARAM_INTERPOLATE_POINT_DATA).GetValue();

        auto filter = vtkSmartPointer<vtkPolyDataPointSampler>::New();
        filter->SetInputData(input_polydata);

        filter->SetDistance(distance);
        filter->SetGenerateVertexPoints(generate_vertex_points);
        filter->SetGenerateEdgePoints(generate_edge_points);
        filter->SetGenerateInteriorPoints(generate_interior_points);
        filter->SetInterpolatePointData(interpolate_point_data);
        filter->SetGenerateVertices(true); // generate cells and not just points, needed for MFX conversion (?)

        filter->Update();

        auto filter_output = filter->GetOutput();
        output_polydata->ShallowCopy(filter_output);
        return kOfxStatOK;
    }
};

// ----------------------------------------------------------------------------

class VtkFeatureEdgesEffect : public VtkEffect {
private:
    const char *PARAM_FEATURE_ANGLE = "FeatureAngle";
    const char *PARAM_FEATURE_EDGES = "FeatureEdges";
    const char *PARAM_BOUNDARY_EDGES = "BoundaryEdges";
    const char *PARAM_NONMANIFOLD_EDGES = "NonManifoldEdges";
    const char *PARAM_MANIFOLD_EDGES = "ManifoldEdges";

public:
    const char* GetName() override {
        return "Feature edges";
    }

    OfxStatus vtkDescribe(OfxParamSetHandle parameters) override {
        AddParam(PARAM_FEATURE_ANGLE, 30).Range(1e-6, 180.0).Label("Feature angle");
        AddParam(PARAM_FEATURE_EDGES, true).Label("Extract feature edges");
        AddParam(PARAM_BOUNDARY_EDGES, false).Label("Extract boundary edges");
        AddParam(PARAM_NONMANIFOLD_EDGES, false).Label("Extract non-manifold edges");
        AddParam(PARAM_MANIFOLD_EDGES, false).Label("Extract manifold edges");
        return kOfxStatOK;
    }

    OfxStatus vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) override {
        double feature_angle = GetParam<double>(PARAM_FEATURE_ANGLE).GetValue();
        bool extract_feature_edges = GetParam<bool>(PARAM_FEATURE_EDGES).GetValue();
        bool extract_boundary_edges = GetParam<bool>(PARAM_BOUNDARY_EDGES).GetValue();
        bool extract_nonmanifold_edges = GetParam<bool>(PARAM_NONMANIFOLD_EDGES).GetValue();
        bool extract_manifold_edges = GetParam<bool>(PARAM_MANIFOLD_EDGES).GetValue();

        auto filter = vtkSmartPointer<vtkFeatureEdges>::New();
        filter->SetInputData(input_polydata);

        filter->SetFeatureAngle(feature_angle);
        filter->SetFeatureEdges(extract_feature_edges);
        filter->SetBoundaryEdges(extract_boundary_edges);
        filter->SetNonManifoldEdges(extract_nonmanifold_edges);
        filter->SetManifoldEdges(extract_manifold_edges);
        filter->SetColoring(false);

        filter->Update();

        // TODO add cleanpolydata to get rid of unused points

        auto filter_output = filter->GetOutput();
        output_polydata->ShallowCopy(filter_output);
        return kOfxStatOK;
    }
};

// ----------------------------------------------------------------------------

class VtkTubeFilterEffect : public VtkEffect {
private:
    const char *PARAM_RADIUS = "Radius";
    const char *PARAM_NUMBER_OF_SIDES = "NumberOfSides";
    const char *PARAM_CAPPING = "Capping";

public:
    const char* GetName() override {
        return "Tube filter";
    }

    OfxStatus vtkDescribe(OfxParamSetHandle parameters) override {
        AddParam(PARAM_RADIUS, 0.05).Range(1e-6, 1e6).Label("Radius");
        AddParam(PARAM_NUMBER_OF_SIDES, 6).Range(3, 1000).Label("Number of sides");
        AddParam(PARAM_CAPPING, true).Label("Cap ends");
        // TODO texture coordinates
        return kOfxStatOK;
    }

    OfxStatus vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) override {
        double radius = GetParam<double>(PARAM_RADIUS).GetValue();
        int number_of_sides = GetParam<int>(PARAM_NUMBER_OF_SIDES).GetValue();
        bool capping = GetParam<bool>(PARAM_CAPPING).GetValue();

        // vtkExtractEdges to create lines even from polygonal mesh
        auto extract_edges_filter = vtkSmartPointer<vtkExtractEdges>::New();
        extract_edges_filter->SetInputData(input_polydata);

        // TODO incorporate optional vtkTubeBender

        // vtkTubeFilter to turn lines into polygonal tubes
        auto tube_filter = vtkSmartPointer<vtkTubeFilter>::New();
        tube_filter->SetInputConnection(extract_edges_filter->GetOutputPort());
        tube_filter->SetRadius(radius);
        tube_filter->SetNumberOfSides(number_of_sides);
        tube_filter->SetCapping(capping);
        tube_filter->SetSidesShareVertices(true);

        // vtkTriangleFilter to convert triangle strips to polygons
        auto triangle_filter = vtkSmartPointer<vtkTriangleFilter>::New();
        triangle_filter->SetInputConnection(tube_filter->GetOutputPort());

        triangle_filter->Update();

        auto filter_output = triangle_filter->GetOutput();
        output_polydata->ShallowCopy(filter_output);
        return kOfxStatOK;
    }
};

// ----------------------------------------------------------------------------

class VtkDecimateProEffect : public VtkEffect {
private:
    const char *PARAM_TARGET_REDUCTION = "TargetReduction";
    const char *PARAM_PRESERVE_TOPOLOGY = "PreserveTopology";
    const char *PARAM_FEATURE_ANGLE = "FeatureAngle";
    const char *PARAM_SPLITTING = "Splitting";
    const char *PARAM_SPLIT_ANGLE = "SplitAngle";
    const char *PARAM_MAXIMUM_ERROR = "MaximumError";
    const char *PARAM_ABSOLUTE_ERROR = "AbsoluteError";
    const char *PARAM_BOUNDARY_VERTEX_DELETION = "BoundaryVertexDeletion";
    const char *PARAM_INFLECTION_POINT_RATIO = "InflectionPointRatio";
    const char *PARAM_DEGREE = "Degree";

public:
    const char* GetName() override {
        return "Decimate (pro)";
    }

    OfxStatus vtkDescribe(OfxParamSetHandle parameters) override {
        AddParam(PARAM_TARGET_REDUCTION, 0.8).Range(1e-6, 1 - 1e-6).Label("Target reduction");
        AddParam(PARAM_PRESERVE_TOPOLOGY, false).Label("Preserve topology");
        AddParam(PARAM_FEATURE_ANGLE, 15.0).Range(0.001, 180.0).Label("Feature angle");
        AddParam(PARAM_SPLITTING, true).Label("Allow splitting");
        AddParam(PARAM_SPLIT_ANGLE, 45.0).Range(0.001, 180.0).Label("Split angle");
        AddParam(PARAM_MAXIMUM_ERROR, 0.01).Range(0, 1e6).Label("Maximum error");
        AddParam(PARAM_ABSOLUTE_ERROR, false).Label("Use absolute error");
        AddParam(PARAM_BOUNDARY_VERTEX_DELETION, true).Label("Allow boundary vertex deletion");
        AddParam(PARAM_INFLECTION_POINT_RATIO, 10.0).Range(1.001, 1e6).Label("Inflection point ratio");
        AddParam(PARAM_DEGREE, 25).Range(3, 1000).Label("Maximum degree of vertex");
        return kOfxStatOK;
    }

    OfxStatus vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) override {
        double target_reduction = GetParam<double>(PARAM_TARGET_REDUCTION).GetValue();
        bool preserve_topology = GetParam<bool>(PARAM_PRESERVE_TOPOLOGY).GetValue();
        double feature_angle = GetParam<double>(PARAM_FEATURE_ANGLE).GetValue();
        bool splitting = GetParam<bool>(PARAM_SPLITTING).GetValue();
        double split_angle = GetParam<double>(PARAM_SPLIT_ANGLE).GetValue();
        double maximum_error = GetParam<double>(PARAM_MAXIMUM_ERROR).GetValue();
        bool absolute_error = GetParam<bool>(PARAM_ABSOLUTE_ERROR).GetValue();
        bool boundary_vertex_deletion = GetParam<bool>(PARAM_BOUNDARY_VERTEX_DELETION).GetValue();
        double inflection_point_ratio = GetParam<double>(PARAM_INFLECTION_POINT_RATIO).GetValue();
        int degree = GetParam<int>(PARAM_DEGREE).GetValue();

        // vtkTriangleFilter to ensure triangle mesh on input
        auto triangle_filter = vtkSmartPointer<vtkTriangleFilter>::New();
        triangle_filter->SetInputData(input_polydata);

        // vtkDecimatePro for main processing
        auto decimate_filter = vtkSmartPointer<vtkDecimatePro>::New();
        decimate_filter->SetInputConnection(triangle_filter->GetOutputPort());
        decimate_filter->SetTargetReduction(target_reduction);
        decimate_filter->SetPreserveTopology(preserve_topology);
        decimate_filter->SetFeatureAngle(feature_angle);
        decimate_filter->SetSplitting(splitting);
        decimate_filter->SetSplitAngle(split_angle);
        decimate_filter->SetMaximumError(maximum_error);
        decimate_filter->SetAbsoluteError(absolute_error);
        decimate_filter->SetBoundaryVertexDeletion(boundary_vertex_deletion);
        decimate_filter->SetInflectionPointRatio(inflection_point_ratio);
        decimate_filter->SetDegree(degree);

        decimate_filter->Update();

        auto filter_output = decimate_filter->GetOutput();
        output_polydata->ShallowCopy(filter_output);
        return kOfxStatOK;
    }
};

// ----------------------------------------------------------------------------

class VtkFillHolesEffect : public VtkEffect {
private:
    const char *PARAM_HOLE_SIZE = "HoleSize";

public:
    const char* GetName() override {
        return "Fill holes";
    }

    OfxStatus vtkDescribe(OfxParamSetHandle parameters) override {
        AddParam(PARAM_HOLE_SIZE, 1.0).Range(1e-6, 1e6).Label("Maximum hole size");
        return kOfxStatOK;
    }

    OfxStatus vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) override {
        double hole_size = GetParam<double>(PARAM_HOLE_SIZE).GetValue();

        auto filter = vtkSmartPointer<vtkFillHolesFilter>::New();
        filter->SetInputData(input_polydata);

        filter->SetHoleSize(hole_size);

        filter->Update();

        auto filter_output = filter->GetOutput();
        output_polydata->ShallowCopy(filter_output);
        return kOfxStatOK;
    }
};

// ----------------------------------------------------------------------------

class VtkIdentityEffect : public VtkEffect {
public:
    const char* GetName() override {
        return "Identity";
    }

    OfxStatus vtkDescribe(OfxParamSetHandle parameters) override {
        return kOfxStatOK;
    }

    OfxStatus vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) override {
        output_polydata->ShallowCopy(input_polydata);
        return kOfxStatOK;
    }
};

// ----------------------------------------------------------------------------

MfxRegister(
        VtkSmoothPolyDataFilterEffect,
        VtkWindowedSincPolyDataFilterEffect,
        VtkPolyDataPointSamplerEffect,
        VtkFeatureEdgesEffect,
        VtkFillHolesEffect,
        VtkTubeFilterEffect,
        VtkDecimateProEffect,
        VtkIdentityEffect
);
