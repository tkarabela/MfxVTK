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

#include "ofxCore.h"
#include "ofxMeshEffect.h"

#include "PluginSupport/MfxRegister"
#include "VtkEffect.h"

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
        AddParam(PARAM_RELAXATION_FACTOR, 0.01).Range(0.0, 1000.0).Label("Relaxation factor");
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
        VtkIdentityEffect
);
