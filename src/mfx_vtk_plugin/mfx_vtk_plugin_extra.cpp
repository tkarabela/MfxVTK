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

#include <vtkCellIterator.h>
#include <vtkPointData.h>
#include <vtkStaticCellLocator.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkFeatureEdges.h>
#include <vtkTriangleFilter.h>
#include <vtkDecimatePro.h>
#include <vtkQuadricClustering.h>
#include <vtkMaskPoints.h>
#include <vtkStaticCellLinks.h>

#include "ofxCore.h"
#include "ofxMeshEffect.h"

#include "PluginSupport/MfxRegister"
#include "VtkEffect.h"
#include "mfx_vtk_utils.h"

// ----------------------------------------------------------------------------

class VtkMaskPointsEffect : public VtkEffect {
private:
    const char *PARAM_RANDOM_MODE = "RandomMode";
    const char *PARAM_RANDOM_MODE_TYPE = "RandomModeType";
    const char *PARAM_ON_RATIO = "OnRatio";
    const char *PARAM_MAXIMUM_NUMBER_OF_POINTS = "MaximumNumberOfPoints";
    // const char *PARAM_RANDOM_SEED = "RandomSeed";

public:
    const char* GetName() override {
        return "Mask points";
    }

    OfxStatus vtkDescribe(OfxParamSetHandle parameters, VtkEffectInputDef &input_mesh, VtkEffectInputDef &output_mesh) override {
        AddParam(PARAM_RANDOM_MODE, true).Label("Use point selection");
        AddParam(PARAM_RANDOM_MODE_TYPE, 0).Range(0, 3).Label("Random distribution type"); // TODO replace this with enum
        AddParam(PARAM_ON_RATIO, 2).Label("Take every n-th point");
        AddParam(PARAM_MAXIMUM_NUMBER_OF_POINTS, 10000).Range(0, 10000000).Label("Maximum number of points");
        // AddParam(PARAM_RANDOM_SEED, 1).Label("Random seed");
        return kOfxStatOK;
    }

    OfxStatus vtkCook(VtkEffectInput &main_input, VtkEffectInput &main_output, std::vector<VtkEffectInput> &extra_inputs) override {
        auto use_random_mode = GetParam<bool>(PARAM_RANDOM_MODE).GetValue();
        auto random_mode_type = GetParam<int>(PARAM_RANDOM_MODE_TYPE).GetValue();
        auto on_ratio = GetParam<int>(PARAM_ON_RATIO).GetValue();
        auto maximum_number_of_points = GetParam<int>(PARAM_MAXIMUM_NUMBER_OF_POINTS).GetValue();
        // auto random_seed = GetParam<int>(PARAM_RANDOM_SEED).GetValue();

        // FIXME until we have eunm, clamp the value manually
        random_mode_type = std::max(0, std::min(3, random_mode_type));

        auto mask_points_filter = vtkSmartPointer<vtkMaskPoints>::New();
        mask_points_filter->SetInputData(main_input.data);

        mask_points_filter->SetRandomMode(use_random_mode);
        mask_points_filter->SetRandomModeType(random_mode_type);
        mask_points_filter->SetOnRatio(on_ratio);
        mask_points_filter->SetMaximumNumberOfPoints(maximum_number_of_points);
        // mask_points_filter->SetRandomSeed(random_seed); // TODO does not exist?

        mask_points_filter->Update();

        auto filter_output = mask_points_filter->GetOutput();
        main_output.data->ShallowCopy(filter_output);
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

    OfxStatus vtkDescribe(OfxParamSetHandle parameters, VtkEffectInputDef &input_mesh, VtkEffectInputDef &output_mesh) override {
        AddParam(PARAM_TARGET_REDUCTION, 0.8).Range(0, 1 - 1e-6).Label("Target reduction");
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

    bool vtkIsIdentity(OfxParamSetHandle parameters) override {
        double target_reduction = GetParam<double>(PARAM_TARGET_REDUCTION).GetValue();
        return !is_positive_double(target_reduction);
    }

    OfxStatus vtkCook(VtkEffectInput &main_input, VtkEffectInput &main_output, std::vector<VtkEffectInput> &extra_inputs) override {
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
        triangle_filter->SetInputData(main_input.data);

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
        main_output.data->ShallowCopy(filter_output);
        return kOfxStatOK;
    }
};

// ----------------------------------------------------------------------------

// TODO replace vtkQuadricClustering with this once it lands post VTK 9.0
//class VtkBinnedDecimationEffect : public VtkEffect {
//private:
//    const char *PARAM_NUMBER_OF_DIVISIONS = "NumberOfDivisions";
//    const char *PARAM_AUTO_ADJUST_NUMBER_OF_DIVISIONS = "AutoAdjustNumberOfDivisions";
//    const char *PARAM_POINT_GENERATION_MODE = "PointGenerationMode";
//
//public:
//    const char* GetName() override {
//        return "Decimate (binned)";
//    }
//
//    OfxStatus vtkDescribe(OfxParamSetHandle parameters) override {
//        AddParam(PARAM_NUMBER_OF_DIVISIONS, std::array<int,3>{256, 256, 256})
//            .Range({2, 2, 2}, {0xffff, 0xffff, 0xffff})
//            .Label("Number of divisions");
//        AddParam(PARAM_AUTO_ADJUST_NUMBER_OF_DIVISIONS, true).Label("Auto adjust number of divisions");
//        AddParam(PARAM_POINT_GENERATION_MODE, 4).Label("Point generation mode (1-4)"); // TODO make this an enum
//        return kOfxStatOK;
//    }
//
//    OfxStatus vtkCook(VtkEffectInput &main_input, VtkEffectInput &main_output, std::vector<VtkEffectInput> &extra_inputs) override {
//        auto number_of_divisions = GetParam<std::array<int,3>>(PARAM_NUMBER_OF_DIVISIONS).GetValue();
//        bool auto_adjust_number_of_divisions = GetParam<bool>(PARAM_AUTO_ADJUST_NUMBER_OF_DIVISIONS).GetValue();
//        bool point_generation_mode = GetParam<int>(PARAM_POINT_GENERATION_MODE).GetValue();
//
//        // TODO do we need this?
//        // vtkTriangleFilter to ensure triangle mesh on input
//        //auto triangle_filter = vtkSmartPointer<vtkTriangleFilter>::New();
//        //triangle_filter->SetInputData(input_polydata);
//
//        // vtkBinnedDecimation for main processing
//        auto decimate_filter = vtkSmartPointer<vtkBinnedDecimation>::New();
//        decimate_filter->SetInputData(input_polydata);
//        decimate_filter->SetNumberOfDivisions(number_of_divisions[0], number_of_divisions[1], number_of_divisions[2]);
//        decimate_filter->SetAutoAdjustNumberOfDivisions(auto_adjust_number_of_divisions);
//        decimate_filter->SetPointGenerationMode(point_generation_mode);
//
//        decimate_filter->Update();
//
//        auto filter_output = decimate_filter->GetOutput();
//        output_polydata->ShallowCopy(filter_output);
//        return kOfxStatOK;
//    }
//};

// ----------------------------------------------------------------------------

class VtkQuadricClusteringEffect : public VtkEffect {
private:
    const char *PARAM_NUMBER_OF_DIVISIONS = "NumberOfDivisions";
    const char *PARAM_AUTO_ADJUST_NUMBER_OF_DIVISIONS = "AutoAdjustNumberOfDivisions";

public:
    const char* GetName() override {
        return "Decimate (quadratic clustering)";
    }

    OfxStatus vtkDescribe(OfxParamSetHandle parameters, VtkEffectInputDef &input_mesh, VtkEffectInputDef &output_mesh) override {
        AddParam(PARAM_NUMBER_OF_DIVISIONS, std::array<int,3>{256, 256, 256})
            .Range({2, 2, 2}, {0xffff, 0xffff, 0xffff})
            .Label("Number of divisions");
        AddParam(PARAM_AUTO_ADJUST_NUMBER_OF_DIVISIONS, true).Label("Auto adjust number of divisions");
        return kOfxStatOK;
    }

    OfxStatus vtkCook(VtkEffectInput &main_input, VtkEffectInput &main_output, std::vector<VtkEffectInput> &extra_inputs) override {
        auto number_of_divisions = GetParam<std::array<int,3>>(PARAM_NUMBER_OF_DIVISIONS).GetValue();
        bool auto_adjust_number_of_divisions = GetParam<bool>(PARAM_AUTO_ADJUST_NUMBER_OF_DIVISIONS).GetValue();

        auto decimate_filter = vtkSmartPointer<vtkQuadricClustering>::New();
        decimate_filter->SetInputData(main_input.data);
        decimate_filter->SetNumberOfDivisions(number_of_divisions[0], number_of_divisions[1], number_of_divisions[2]);
        decimate_filter->SetAutoAdjustNumberOfDivisions(auto_adjust_number_of_divisions);

        decimate_filter->Update();

        auto filter_output = decimate_filter->GetOutput();
        main_output.data->ShallowCopy(filter_output);
        return kOfxStatOK;
    }
};

// ----------------------------------------------------------------------------

class VtkIdentityEffect : public VtkEffect {
public:
    const char *PARAM_ACTION_IS_IDENTITY = "ActionIsIdentity";

    const char* GetName() override {
        return "Identity";
    }

    OfxStatus vtkDescribe(OfxParamSetHandle parameters, VtkEffectInputDef &input_mesh, VtkEffectInputDef &output_mesh) override {
        AddParam(PARAM_ACTION_IS_IDENTITY, false).Label("kOfxMeshEffectActionIsIdentity");
        return kOfxStatOK;
    }

    bool vtkIsIdentity(OfxParamSetHandle parameters) override {
        bool action_is_identity = GetParam<bool>(PARAM_ACTION_IS_IDENTITY).GetValue();
        return action_is_identity;
    }

    OfxStatus vtkCook(VtkEffectInput &main_input, VtkEffectInput &main_output, std::vector<VtkEffectInput> &extra_inputs) override {
        main_output.data->ShallowCopy(main_input.data);
        return kOfxStatOK;
    }
};

// ----------------------------------------------------------------------------

static void print_matrix(double m[16]) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            printf("\t%.4g", m[4*i+j]);
        }
        printf("\n");
    }
}

static void copy_matrix(double* in, double *out) {
    for (int i = 0; i < 16; ++i) {
        out[i] = in[i];
    }
}

static void invert_matrix(double m[16]) {
    double **A = new double*[4];
    double **AI = new double*[4];
    for (int i = 0; i < 4; ++i) {
        A[i] = new double[4];
        AI[i] = new double[4];
        for (int j = 0; j < 4; ++j) {
            A[i][j] = m[4*i + j];
            AI[i][j] = 0.0;
        }
    }

    vtkMath::InvertMatrix(A, AI, 4);

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            m[4*i + j] = AI[i][j];
        }
    }

    for (int i = 0; i < 4; ++i) {
        delete[] A[i];
        delete[] AI[i];
    }

    delete[] A;
    delete[] AI;
}

static void multiply_matrix(double A[16], double B[16], double *C) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            C[4*i+j] = 0.0;
            for (int k = 0; k < 4; ++k) {
                C[4*i+j] += A[4*i+k] * B[4*k+j];
            }
        }
    }
}

static void transform_coordinates(double A[16], double xyz[3]) {
    double x[4] = { xyz[0], xyz[1], xyz[2], 1.0 }, b[4] = {0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            b[i] += A[4*i+j] * x[j];
        }
    }
    for (int i = 0; i < 3; ++i) {
        xyz[i] = b[i];
    }
}

#if 0
class AppendInputEffect : public MfxEffect {
public:
    const char *SECOND_INPUT_NAME = "SecondInput";

    const char* GetName() override {
        return "Append input (second mesh)";
    }

    OfxStatus Describe(OfxMeshEffectHandle descriptor) override {
        auto main_input = AddInput(kOfxMeshMainInput);
        main_input.RequestTransform(true);
        AddInput(kOfxMeshMainOutput);
        auto second_mesh = AddInput(SECOND_INPUT_NAME);
        second_mesh.Label("Second input");
        second_mesh.RequestTransform(true);
        return kOfxStatOK;
    }

    OfxStatus Cook(OfxMeshEffectHandle instance) override {
        MfxMesh input_mesh = GetInput(kOfxMeshMainInput).GetMesh();
        MfxMesh second_mesh = GetInput(SECOND_INPUT_NAME).GetMesh();
        MfxMesh output_mesh = GetInput(kOfxMeshMainOutput).GetMesh();

        MfxMesh* input_meshes[2] = { &input_mesh, &second_mesh };

        int pointCount = 0, vertexCount = 0, faceCount = 0, noLooseEdge = 1, constantFaceCount = -1;

        for (int i = 0; i < 2; i++) {
            if (!input_meshes[i]->IsValid()) {
                printf("Error - input mesh %d is not valid\n", i);
                return kOfxStatFailed; // XXX ??
            }

            MfxMeshProps props;
            input_meshes[i]->FetchProperties(props);

            printf("OFX input_mesh[%d] has %d points %d vertices %d faces\n", i, props.pointCount, props.vertexCount, props.faceCount);

            pointCount += props.pointCount;
            vertexCount += props.vertexCount;
            faceCount += props.faceCount;
            if (!props.noLooseEdge) noLooseEdge = 0;
        }

        double *m1 = nullptr;
        double *m2 = nullptr;
        input_mesh.FetchTransform(reinterpret_cast<double **>(&m1));
        second_mesh.FetchTransform(reinterpret_cast<double **>(&m2));

        printf("m1:\n");
        print_matrix(m1);
        printf("\n");
        printf("m2:\n");
        print_matrix(m2);
        printf("\n");

        double m1_inv[16];
        copy_matrix(m1, m1_inv);
        invert_matrix(m1_inv);

        printf("inv(m1):\n");
        print_matrix(m1_inv);
        printf("\n");

        double m1_inv_m2[16];
        multiply_matrix(m1_inv, m2, m1_inv_m2);
        printf("inv(m1) . m2:\n");
        print_matrix(m1_inv_m2);
        printf("\n");

        printf("OFX output mesh has %d points %d vertices %d faces\n", pointCount, vertexCount, faceCount);
        output_mesh.Allocate(pointCount,
                             vertexCount,
                             faceCount,
                             noLooseEdge,
                             constantFaceCount);

        MfxMeshProps output_mesh_properties;
        output_mesh.FetchProperties(output_mesh_properties);

        MfxAttribute output_pointPos = output_mesh.GetPointAttribute(kOfxMeshAttribPointPosition);
        MfxAttribute output_vertPoint = output_mesh.GetVertexAttribute(kOfxMeshAttribVertexPoint);
        MfxAttribute output_faceLen = output_mesh.GetFaceAttribute(kOfxMeshAttribFaceCounts);

        int i_pointPos = 0, i_vertPoint = 0, i_faceLen = 0, point_start_id = 0;

        for (int i = 0; i < 2; i++) {
            MfxMeshProps input_mesh_props;
            input_meshes[i]->FetchProperties(input_mesh_props);
            int constant_face_count = input_mesh_props.constantFaceCount;

            MfxAttribute input_pointPos = input_meshes[i]->GetPointAttribute(kOfxMeshAttribPointPosition);
            MfxAttribute input_vertPoint = input_meshes[i]->GetVertexAttribute(kOfxMeshAttribVertexPoint);
            MfxAttribute input_faceLen = input_meshes[i]->GetFaceAttribute(kOfxMeshAttribFaceCounts);

            MfxAttributeProps in_prop, out_prop;

            input_pointPos.FetchProperties(in_prop);
            output_pointPos.FetchProperties(out_prop);
            for (int j = 0; j < input_mesh_props.pointCount; j++) {
                double xyz[3];
                for (int k = 0; k < 3; k++) {
                    float x;
                    in_prop.GetValue(&x, j, k);
                    xyz[k] = x;
                    // out_prop.SetValue(x, i_pointPos, k);
                }

                if (i > 0) {
                    transform_coordinates(m1_inv_m2, xyz);
                }

                for (int k = 0; k < 3; k++) {
                    float x = xyz[k];
                    out_prop.SetValue(x, i_pointPos, k);
                }

                i_pointPos++;
            }

            input_vertPoint.FetchProperties(in_prop);
            output_vertPoint.FetchProperties(out_prop);
            for (int j = 0; j < input_mesh_props.vertexCount; j++) {
                int x;
                in_prop.GetValue(&x, j);
                //out_prop.SetValue(x, i_vertPoint);
                out_prop.SetValue(x + point_start_id, i_vertPoint);
                //printf("output_vertPoint[%d] = %d\n", i_vertPoint, x+vertex_start_id);
                i_vertPoint++;
            }
            point_start_id += input_mesh_props.pointCount;

            input_faceLen.FetchProperties(in_prop);
            output_faceLen.FetchProperties(out_prop);
            for (int j = 0; j < input_mesh_props.faceCount; j++) {
                int x;
                if (constant_face_count > 0) {
                    x = constant_face_count;
                } else {
                    in_prop.GetValue(&x, j);
                }
                out_prop.SetValue(x, i_faceLen);
                i_faceLen++;
            }
        }

        output_mesh.Release();
        input_mesh.Release();
        second_mesh.Release();

        return kOfxStatOK;
    }
};
#endif

// ----------------------------------------------------------------------------

class IdentityForwardAttributesEffect : public MfxEffect {
public:
    const char* GetName() override {
        return "Identity (forward attributes)";
    }

    OfxStatus Describe(OfxMeshEffectHandle descriptor) override {
        AddInput(kOfxMeshMainInput);
        AddInput(kOfxMeshMainOutput);
        return kOfxStatOK;
    }

    OfxStatus Cook(OfxMeshEffectHandle instance) override {
        MfxMesh input_mesh = GetInput(kOfxMeshMainInput).GetMesh();
        MfxMesh output_mesh = GetInput(kOfxMeshMainOutput).GetMesh();

        MfxMeshProps input_mesh_properties, output_mesh_properties;
        input_mesh.FetchProperties(input_mesh_properties);
        output_mesh.FetchProperties(output_mesh_properties);

        MfxAttribute input_pointPos = input_mesh.GetPointAttribute(kOfxMeshAttribPointPosition);
        MfxAttribute input_vertPoint = input_mesh.GetVertexAttribute(kOfxMeshAttribVertexPoint);
        MfxAttribute input_faceLen = input_mesh.GetFaceAttribute(kOfxMeshAttribFaceCounts);

        MfxAttribute output_pointPos = output_mesh.GetPointAttribute(kOfxMeshAttribPointPosition);
        MfxAttribute output_vertPoint = output_mesh.GetVertexAttribute(kOfxMeshAttribVertexPoint);
        MfxAttribute output_faceLen = output_mesh.GetFaceAttribute(kOfxMeshAttribFaceCounts);

        output_pointPos.ForwardFrom(input_pointPos);
        output_vertPoint.ForwardFrom(input_vertPoint);
        output_faceLen.ForwardFrom(input_faceLen);

        output_mesh.Allocate(input_mesh_properties.pointCount,
                             input_mesh_properties.vertexCount,
                             input_mesh_properties.faceCount,
                             input_mesh_properties.noLooseEdge,
                             input_mesh_properties.constantFaceCount);

        output_mesh.Release();
        input_mesh.Release();

        return kOfxStatOK;
    }
};

// ----------------------------------------------------------------------------

MfxRegister(
    // these are unsorted effects wrapped directly from VTK
    VtkMaskPointsEffect,
    VtkDecimateProEffect,
    VtkQuadricClusteringEffect,

    // these effects are interesting only for development of Open Mesh Effect
    VtkIdentityEffect,
    IdentityForwardAttributesEffect
    //AppendInputEffect
);
