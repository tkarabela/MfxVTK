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
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>
#include <vtkFeatureEdges.h>
#include <vtkTriangleFilter.h>
#include <vtkDecimatePro.h>
#include <vtkQuadricClustering.h>
#include <vtkMaskPoints.h>
#include <vtkStaticCellLinks.h>
#include <vtkThreshold.h>
#include <vtkGeometryFilter.h>
#include <vtkIdFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkModifiedBSPTree.h>
#include <chrono>
#include <vtkOctreePointLocator.h>

#include "ofxCore.h"
#include "ofxMeshEffect.h"

#include "PluginSupport/MfxRegister"
#include "VtkEffect.h"
#include "mfx_vtk_utils.h"

// ----------------------------------------------------------------------------

class VtkPokeEffect : public VtkEffect {
private:
    const char *PARAM_MAX_DISTANCE = "MaxDistance";
    const char *PARAM_FALLOFF_RADIUS = "FalloffRadius";
    const char *PARAM_FALLOFF_EXPONENT = "FalloffExponent";
    const char *PARAM_OFFSET = "Offset";

public:
    const char* GetName() override {
        return "Poke";
    }

    OfxStatus vtkDescribe(OfxParamSetHandle parameters) override {
        // TODO support collision with other mesh
        // AddInput("Collider");

        // TODO support specifying attribute for using main mesh as collider

        AddParam(PARAM_MAX_DISTANCE, 0.0).Range(0.0, 1e6).Label("Maximum distance");
        AddParam(PARAM_FALLOFF_RADIUS, 0.01).Range(0.0, 1e6).Label("Falloff radius");
        AddParam(PARAM_FALLOFF_EXPONENT, 2.0).Range(1e-2, 1e2).Label("Falloff exponent");
        AddParam(PARAM_OFFSET, 0.0).Range(0.0, 1e6).Label("Offset");
        return kOfxStatOK;
    }

    OfxStatus vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) override {
        auto max_distance = GetParam<double>(PARAM_MAX_DISTANCE).GetValue();
        auto falloff_radius = GetParam<double>(PARAM_FALLOFF_RADIUS).GetValue();
        auto falloff_exponent = GetParam<double>(PARAM_FALLOFF_EXPONENT).GetValue();
        auto offset = GetParam<double>(PARAM_OFFSET).GetValue();

        if (!input_polydata->GetPointData()->HasArray("color0")) {
            printf("MfxVTK - VtkPokeEffect - error, input must have vertex color attribute called 'color0'\n");
            return kOfxStatErrValue;
        }

        // XXX NOTE TO USER - paint mesh white and collider black...
        // TODO we really want to support gray colors as well, to exclude boundary of collider and mesh


        auto t0 = std::chrono::system_clock::now();

        // remember point IDs so that we can map deformed mesh onto input_polydata
        const char *point_id_array_name = "_PointId";
        auto id_filter = vtkSmartPointer<vtkIdFilter>::New();
        id_filter->SetInputData(input_polydata);
        id_filter->SetPointIds(true);
        id_filter->SetPointIdsArrayName(point_id_array_name);
        id_filter->SetCellIds(false);

        // get normals
        // TODO read normals from input instead of recalculating (but make sure we have cell normals in collider)
        auto normals_filter = vtkSmartPointer<vtkPolyDataNormals>::New();
        normals_filter->SetInputConnection(id_filter->GetOutputPort());
        normals_filter->SetSplitting(false);
        normals_filter->SetNonManifoldTraversal(true);
        normals_filter->SetComputePointNormals(true);
        normals_filter->SetComputeCellNormals(true); // for collider
        normals_filter->SetAutoOrientNormals(true);

        // separate input into collider part and mesh part
        auto threshold_filter = vtkSmartPointer<vtkThreshold>::New();
        threshold_filter->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "color0");
        threshold_filter->ThresholdByUpper(0.5);
        threshold_filter->SetInvert(true);
        threshold_filter->SetInputConnection(normals_filter->GetOutputPort());

        // vtkThreshold outputs unstructured grid, get vtkPolyData again
        auto geometry_filter = vtkSmartPointer<vtkGeometryFilter>::New();
        geometry_filter->SetInputConnection(threshold_filter->GetOutputPort());
        geometry_filter->SetMerging(false);
        geometry_filter->Update();

        auto mesh_polydata = vtkSmartPointer<vtkPolyData>::New();
        mesh_polydata->ShallowCopy(geometry_filter->GetOutput());

        threshold_filter->SetInvert(false);
        threshold_filter->Update();
        geometry_filter->Update();
        auto collider_polydata = vtkSmartPointer<vtkPolyData>::New();
        collider_polydata->ShallowCopy(geometry_filter->GetOutput());

        auto t1 = std::chrono::system_clock::now();

        // do the deformation
        // TODO crop mesh_polydata by collider bounding geometry, to reduce cost?
        vtkCookInner(mesh_polydata, collider_polydata, max_distance, falloff_radius, falloff_exponent, offset);

        auto t2 = std::chrono::system_clock::now();

        // map deformation back onto main mesh
        output_polydata->ShallowCopy(input_polydata);
        output_polydata->GetPoints()->DeepCopy(input_polydata->GetPoints()); // maybe skip if we're okay changing input in place

        auto mesh_id_array = mesh_polydata->GetPointData()->GetArray(point_id_array_name);
        auto mesh_points = mesh_polydata->GetPoints();
        auto output_points = output_polydata->GetPoints();
        // TODO parallelize this
        for (int pid_mesh = 0; pid_mesh < mesh_polydata->GetNumberOfPoints(); pid_mesh++) {
            int pid_output = mesh_id_array->GetVariantValue(pid_mesh).ToInt();
            double p[3];
            mesh_points->GetPoint(pid_mesh, p);
            output_points->SetPoint(pid_output, p);
        }

        auto t3 = std::chrono::system_clock::now();

        auto dt = [](std::chrono::system_clock::time_point t1, std::chrono::system_clock::time_point t2) -> int {
            return std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        };

        printf("\n\tVtkPokeEffect timings [ms] dt1=%d, dt2=%d, dt3=%d\n",
               dt(t0, t1), dt(t1, t2), dt(t2, t3));

        return kOfxStatOK;
    }

    static void vtkCookInner(vtkPolyData *mesh_polydata, vtkPolyData *collider_polydata,
                             double max_distance, double falloff_radius, double falloff_exponent, double offset) {
        if (!is_positive_double(max_distance)) {
            max_distance = mesh_polydata->GetLength(); // bounding box diagonal
        }

        auto mesh_normals = mesh_polydata->GetPointData()->GetArray("Normals");
        auto collider_cell_normals = collider_polydata->GetCellData()->GetArray("Normals");

        auto new_mesh_points = vtkSmartPointer<vtkPoints>::New();
        new_mesh_points->DeepCopy(mesh_polydata->GetPoints()); // copy this so iterations do not affect each other

        auto displacement_distance_array = vtkSmartPointer<vtkFloatArray>::New();
        displacement_distance_array->SetNumberOfComponents(1);
        displacement_distance_array->SetNumberOfTuples(mesh_polydata->GetNumberOfPoints());
        displacement_distance_array->FillValue(0.0);

        auto mesh_cell_locator = vtkSmartPointer<vtkModifiedBSPTree>::New();
        mesh_cell_locator->SetDataSet(mesh_polydata);
        mesh_cell_locator->Update();

        auto collider_cell_locator = vtkSmartPointer<vtkModifiedBSPTree>::New();
        collider_cell_locator->SetDataSet(collider_polydata);
        collider_cell_locator->Update();

        auto tmp_points_mesh = vtkSmartPointer<vtkPoints>::New();
        auto tmp_points_collider = vtkSmartPointer<vtkPoints>::New();
        auto tmp_cells_collider = vtkSmartPointer<vtkIdList>::New();

        // bounding box of affected mesh area
        double min_coord[3] = { DBL_MAX, DBL_MAX, DBL_MAX };
        double max_coord[3] = { DBL_MIN, DBL_MIN, DBL_MIN };

        double collider_bounds[6];
        collider_polydata->GetBounds(collider_bounds);

        double min_coord_estimate[3] = { collider_bounds[0] - max_distance,
                                         collider_bounds[2] - max_distance,
                                         collider_bounds[4] - max_distance };
        double max_coord_estimate[3] = { collider_bounds[1] + max_distance,
                                         collider_bounds[2] + max_distance,
                                         collider_bounds[5] + max_distance };

        // step 1 - clear mesh collision with collider
        for (int i = 0; i < mesh_polydata->GetNumberOfPoints(); i++) {
            double p[3];
            double mesh_point[3], collider_point[3];
            double mesh_normal[3], collider_cell_normal[3];
            mesh_polydata->GetPoint(i, p);

            // early exit, mesh point too far away from collider
            if (p[0] < min_coord_estimate[0] || p[0] > max_coord_estimate[0] ||
                p[1] < min_coord_estimate[1] || p[1] > max_coord_estimate[1] ||
                p[2] < min_coord_estimate[2] || p[2] > max_coord_estimate[2]) {
                break;
            }

            mesh_normals->GetTuple(i, mesh_normal);

            const double p0[3] =    { p[0] - 1e-6*mesh_normal[0],
                                      p[1] - 1e-6*mesh_normal[1],
                                      p[2] - 1e-6*mesh_normal[2] }; // probe line start, just under point
            const double p_inf[3] = { p[0] - max_distance*mesh_normal[0],
                                      p[1] - max_distance*mesh_normal[1],
                                      p[2] - max_distance*mesh_normal[2] }; // probe line end, deep under point
            const double tol = 1e-4; // ???

            tmp_points_collider->Reset();
            tmp_cells_collider->Reset();
            collider_cell_locator->IntersectWithLine(p0, p_inf, tol, tmp_points_collider, tmp_cells_collider);
            if (tmp_points_collider->GetNumberOfPoints() == 0) {
                continue; // point does not have collider "behind" it
            }
            tmp_points_mesh->Reset();
            mesh_cell_locator->IntersectWithLine(p0, p_inf, tol, tmp_points_mesh, nullptr);

            tmp_points_collider->GetPoint(0, collider_point);
            double this_to_collider_sq = vec3_squared_distance(p, collider_point);
            double this_to_mesh_sq = DBL_MAX;
            double mesh_to_collider_sq = DBL_MAX;
            if (tmp_points_mesh->GetNumberOfPoints() > 0) {
                tmp_points_mesh->GetPoint(0, mesh_point);
                this_to_mesh_sq = vec3_squared_distance(p, mesh_point);

                // mesh_to_collider_sq := distance from first mesh intersection to preceding collider intersection
                for (int j = 0; j < tmp_points_collider->GetNumberOfPoints(); j++) {
                    tmp_points_collider->GetPoint(j, collider_point);
                    double tmp_collider_distance_sq = vec3_squared_distance(p, collider_point);
                    if (tmp_collider_distance_sq > this_to_mesh_sq) break;
                    mesh_to_collider_sq = std::abs(tmp_collider_distance_sq - this_to_mesh_sq);
                }
            }

            // step 1, final breakdown
            if (this_to_mesh_sq < this_to_collider_sq) {
                // collider is outside this section of mesh
                //printf("point %d - NO, collider outside\n", i);
            } else if (this_to_collider_sq > mesh_to_collider_sq) {
                // other side of mesh is closer to the collider
                //printf("point %d - NO, other side closer to collider\n", i);
            } else {
                collider_cell_normals->GetTuple(tmp_cells_collider->GetId(0), collider_cell_normal);
                double x = vec3_dot(mesh_normal, collider_cell_normal);
                
                if (x > 0) {
                    // bad angle, collider and mesh point look the same way
                    //printf("point %d - NO, bad angle %g\n", i, x);
                } else {
                    double distance = std::sqrt(this_to_collider_sq);
                    /*
                    double new_p[3] = { p[0] - mesh_normal[0]*distance,
                                        p[1] - mesh_normal[1]*distance,
                                        p[2] - mesh_normal[2]*distance };
                    new_mesh_points->SetPoint(i, new_p);
                    */
                    displacement_distance_array->SetValue(i, distance);
                    vec3_min(min_coord, p);
                    vec3_max(max_coord, p);
                    //printf("point %d - OK, moving by distance %g\n", i, distance);
                }
            }
        }

        // step 2 - handle falloff
        auto mesh_falloff_locator = vtkSmartPointer<vtkOctreePointLocator>::New();
        mesh_falloff_locator->SetDataSet(mesh_polydata);
        mesh_falloff_locator->Update();

        auto points_in_area_of_effect = vtkSmartPointer<vtkIdTypeArray>::New();
        double area_param[6] = { min_coord[0] - falloff_radius, max_coord[0] + falloff_radius,
                                 min_coord[1] - falloff_radius, max_coord[1] + falloff_radius,
                                 min_coord[2] - falloff_radius, max_coord[2] + falloff_radius };

        mesh_falloff_locator->FindPointsInArea(area_param, points_in_area_of_effect);
        auto tmp_points_falloff = vtkSmartPointer<vtkIdList>::New();

        const double falloff_radius_sq = falloff_radius*falloff_radius;

        for (int _i = 0; _i < points_in_area_of_effect->GetNumberOfValues(); _i++) {
            int i = points_in_area_of_effect->GetValue(_i);
            double p0[3], mesh_normal[3];
            mesh_polydata->GetPoint(i, p0);
            mesh_normals->GetTuple(i, mesh_normal);

            tmp_points_falloff->Reset();

            mesh_falloff_locator->FindPointsWithinRadius(falloff_radius, p0, tmp_points_falloff);
            int num_neighbors = tmp_points_falloff->GetNumberOfIds();
            double weight_sum = 0.0;
            double weighted_displacement_sum = 0.0;

            // XXX we should really traverse the manifold instead of doing ball search like this
            // XXX we should really average normals after deformation or sth instead of using previous normals
            // XXX can this be parallel?
            for (int _j = 0; _j < num_neighbors; _j++) {
                int j = tmp_points_falloff->GetId(_j);
                double p[3];
                mesh_polydata->GetPoint(j, p);
                double displacement = displacement_distance_array->GetValue(j);
                double distance_sq = vec3_squared_distance(p, p0);
                double radius_ratio = distance_sq / falloff_radius_sq;
                double weight = std::pow(radius_ratio, falloff_exponent);

                weight_sum += weight;
                weighted_displacement_sum += /* weight * */ displacement; // XXX incorrect but looks better
            }

            double weighted_average_displacement_distance = weighted_displacement_sum / weight_sum;
            double old_distance = displacement_distance_array->GetValue(i);
            double new_distance = std::min(max_distance, std::max(old_distance, weighted_average_displacement_distance));

            double new_p[3] = { p0[0] - mesh_normal[0]*new_distance,
                                p0[1] - mesh_normal[1]*new_distance,
                                p0[2] - mesh_normal[2]*new_distance };
            new_mesh_points->SetPoint(i, new_p);
        }

        // final step - overwrite original coordinates with new ones
        mesh_polydata->GetPoints()->ShallowCopy(new_mesh_points);
    }
};

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

    OfxStatus vtkDescribe(OfxParamSetHandle parameters) override {
        AddParam(PARAM_RANDOM_MODE, true).Label("Use point selection");
        AddParam(PARAM_RANDOM_MODE_TYPE, 0).Range(0, 3).Label("Random distribution type"); // TODO replace this with enum
        AddParam(PARAM_ON_RATIO, 2).Label("Take every n-th point");
        AddParam(PARAM_MAXIMUM_NUMBER_OF_POINTS, 10000).Range(0, 10000000).Label("Maximum number of points");
        // AddParam(PARAM_RANDOM_SEED, 1).Label("Random seed");
        return kOfxStatOK;
    }

    OfxStatus vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) override {
        auto use_random_mode = GetParam<bool>(PARAM_RANDOM_MODE).GetValue();
        auto random_mode_type = GetParam<int>(PARAM_RANDOM_MODE_TYPE).GetValue();
        auto on_ratio = GetParam<int>(PARAM_ON_RATIO).GetValue();
        auto maximum_number_of_points = GetParam<int>(PARAM_MAXIMUM_NUMBER_OF_POINTS).GetValue();
        // auto random_seed = GetParam<int>(PARAM_RANDOM_SEED).GetValue();

        // FIXME until we have eunm, clamp the value manually
        random_mode_type = std::max(0, std::min(3, random_mode_type));

        auto mask_points_filter = vtkSmartPointer<vtkMaskPoints>::New();
        mask_points_filter->SetInputData(input_polydata);

        mask_points_filter->SetRandomMode(use_random_mode);
        mask_points_filter->SetRandomModeType(random_mode_type);
        mask_points_filter->SetOnRatio(on_ratio);
        mask_points_filter->SetMaximumNumberOfPoints(maximum_number_of_points);
        // mask_points_filter->SetRandomSeed(random_seed); // TODO does not exist?

        mask_points_filter->Update();

        auto filter_output = mask_points_filter->GetOutput();
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
//    OfxStatus vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) override {
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

    OfxStatus vtkDescribe(OfxParamSetHandle parameters) override {
        AddParam(PARAM_NUMBER_OF_DIVISIONS, std::array<int,3>{256, 256, 256})
            .Range({2, 2, 2}, {0xffff, 0xffff, 0xffff})
            .Label("Number of divisions");
        AddParam(PARAM_AUTO_ADJUST_NUMBER_OF_DIVISIONS, true).Label("Auto adjust number of divisions");
        return kOfxStatOK;
    }

    OfxStatus vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) override {
        auto number_of_divisions = GetParam<std::array<int,3>>(PARAM_NUMBER_OF_DIVISIONS).GetValue();
        bool auto_adjust_number_of_divisions = GetParam<bool>(PARAM_AUTO_ADJUST_NUMBER_OF_DIVISIONS).GetValue();

        auto decimate_filter = vtkSmartPointer<vtkQuadricClustering>::New();
        decimate_filter->SetInputData(input_polydata);
        decimate_filter->SetNumberOfDivisions(number_of_divisions[0], number_of_divisions[1], number_of_divisions[2]);
        decimate_filter->SetAutoAdjustNumberOfDivisions(auto_adjust_number_of_divisions);

        decimate_filter->Update();

        auto filter_output = decimate_filter->GetOutput();
        output_polydata->ShallowCopy(filter_output);
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

    OfxStatus vtkDescribe(OfxParamSetHandle parameters) override {
        AddParam(PARAM_ACTION_IS_IDENTITY, false).Label("kOfxMeshEffectActionIsIdentity");
        return kOfxStatOK;
    }

    bool vtkIsIdentity(OfxParamSetHandle parameters) override {
        bool action_is_identity = GetParam<bool>(PARAM_ACTION_IS_IDENTITY).GetValue();
        return action_is_identity;
    }

    OfxStatus vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) override {
        output_polydata->ShallowCopy(input_polydata);
        return kOfxStatOK;
    }
};

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
    // effects in development
    VtkPokeEffect,

    // these are unsorted effects wrapped directly from VTK
    VtkMaskPointsEffect,
    VtkDecimateProEffect,
    VtkQuadricClusteringEffect,

    // these effects are interesting only for development of Open Mesh Effect
    VtkIdentityEffect,
    IdentityForwardAttributesEffect
);
