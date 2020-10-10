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
        AddParam(PARAM_CONVERGENCE, 0.0).Range(0.0, 1000.0).Label("Convergence");
        AddParam(PARAM_BOUNDARY_SMOOTHING, true).Label("Boundary smoothing");
        AddParam(PARAM_RELAXATION_FACTOR, 0.01).Range(0.0, 1000.0).Label("Relaxation factor");
        AddParam(PARAM_FEATURE_EDGE_SMOOTHING, false).Label("Feature edge smoothing");
        AddParam(PARAM_FEATURE_ANGLE, 45.0).Range(0.001, 180.0).Label("Feature angle");
        AddParam(PARAM_EDGE_ANGLE, 15.0).Range(0.001, 180.0).Label("Edge angle");
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
        VtkIdentityEffect
);
