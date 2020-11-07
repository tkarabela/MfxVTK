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

#include "VtkPokeEffect.h"
#include "mfx_vtk_utils.h"
#include <chrono>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkIdFilter.h>
#include <vtkIdTypeArray.h>
#include <vtkFloatArray.h>
#include <vtkPolyDataNormals.h>
#include <vtkThreshold.h>
#include <vtkGeometryFilter.h>
#include <vtkModifiedBSPTree.h>
#include <vtkOctreePointLocator.h>

const char *VtkPokeEffect::GetName() {
    return "Poke";
}

OfxStatus VtkPokeEffect::vtkDescribe(OfxParamSetHandle parameters) {
    // TODO support collision with other mesh
    // AddInput("Collider");

    // TODO support specifying attribute for using main mesh as collider

    AddParam(PARAM_MAX_DISTANCE, 0.0).Range(0.0, 1e6).Label("Maximum distance");
    AddParam(PARAM_FALLOFF_RADIUS, 0.01).Range(0.0, 1e6).Label("Falloff radius");
    AddParam(PARAM_FALLOFF_EXPONENT, 2.0).Range(1e-2, 1e2).Label("Falloff exponent");
    AddParam(PARAM_OFFSET, 0.0).Range(0.0, 1e6).Label("Offset");
    AddParam(PARAM_DEBUG, false).Label("Debug");
    return kOfxStatOK;
}

OfxStatus VtkPokeEffect::vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) {
    auto max_distance = GetParam<double>(PARAM_MAX_DISTANCE).GetValue();
    auto falloff_radius = GetParam<double>(PARAM_FALLOFF_RADIUS).GetValue();
    auto falloff_exponent = GetParam<double>(PARAM_FALLOFF_EXPONENT).GetValue();
    auto offset = GetParam<double>(PARAM_OFFSET).GetValue();
    auto debug = GetParam<bool>(PARAM_DEBUG).GetValue();

    if (!input_polydata->GetPointData()->HasArray("color0")) {
        printf("VtkPokeEffect - error, input must have vertex color attribute called 'color0'\n");
        return kOfxStatFailed;
    }

    // XXX NOTE TO USER - paint mesh white and collider black...

    vtkCook_inner(input_polydata, output_polydata, max_distance, falloff_radius, falloff_exponent, offset, debug);

    return kOfxStatOK;
}

OfxStatus VtkPokeEffect::vtkCook_inner(vtkPolyData *input_polydata, vtkPolyData *output_polydata, double max_distance,
                                       double falloff_radius, double falloff_exponent, double offset, bool debug) {
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
    threshold_filter->ThresholdBetween(255, 255);
    threshold_filter->SetInputConnection(normals_filter->GetOutputPort());

    // vtkThreshold outputs unstructured grid, get vtkPolyData again
    auto geometry_filter = vtkSmartPointer<vtkGeometryFilter>::New();
    geometry_filter->SetInputConnection(threshold_filter->GetOutputPort());
    geometry_filter->SetMerging(false);
    geometry_filter->Update();

    auto mesh_polydata = vtkSmartPointer<vtkPolyData>::New();
    mesh_polydata->ShallowCopy(geometry_filter->GetOutput());

    threshold_filter->ThresholdBetween(0, 0);
    threshold_filter->Update();
    geometry_filter->Update();
    auto collider_polydata = vtkSmartPointer<vtkPolyData>::New();
    collider_polydata->ShallowCopy(geometry_filter->GetOutput());

    auto t1 = std::chrono::system_clock::now();

    // do the deformation
    printf("\n\ttotal elements: %d  collider: %d   mesh: %d\n", input_polydata->GetNumberOfCells(), collider_polydata->GetNumberOfCells(), mesh_polydata->GetNumberOfCells());

    // TODO crop mesh_polydata by collider bounding geometry, to reduce cost?
    if (!is_positive_double(max_distance)) {
        max_distance = mesh_polydata->GetLength(); // bounding box diagonal
    }
    printf("VtkPokeEffect - max_distance = %g\n", max_distance);

    auto contacts = evaluate_collision(mesh_polydata, collider_polydata, max_distance, offset);
    auto t2 = std::chrono::system_clock::now();
    handle_reaction_simple(mesh_polydata, falloff_radius, falloff_exponent, contacts, debug, max_distance);
    auto t3 = std::chrono::system_clock::now();


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

    auto t4 = std::chrono::system_clock::now();

    auto dt = [](std::chrono::system_clock::time_point t1, std::chrono::system_clock::time_point t2) -> int {
        return std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    };

    printf("\n\tVtkPokeEffect timings [ms] dt1=%d, dt2(collision)=%d, dt3(reaction)=%d, dt4=%d\n",
           dt(t0, t1), dt(t1, t2), dt(t2, t3), dt(t3, t4));

    return kOfxStatOK;
}

std::vector<VtkPokeEffect::Contact>
VtkPokeEffect::evaluate_collision(vtkPolyData *mesh_polydata, vtkPolyData *collider_polydata, double max_distance,
                                  double offset) {
    std::vector<Contact> contacts;

    // FIXME max_distance has some trouble when its not set to large value

    if (offset > 0) {
        // TODO handle offset
        printf("VtkPokeEffect - warning, nonzero offset not handled yet\n");
    }

    auto mesh_normals = mesh_polydata->GetPointData()->GetArray("Normals");
    auto collider_cell_normals = collider_polydata->GetCellData()->GetArray("Normals");

    auto new_mesh_points = vtkSmartPointer<vtkPoints>::New();
    new_mesh_points->DeepCopy(mesh_polydata->GetPoints()); // copy this so iterations do not affect each other

    auto mesh_cell_locator = vtkSmartPointer<vtkModifiedBSPTree>::New();
    mesh_cell_locator->SetDataSet(mesh_polydata);
    mesh_cell_locator->Update();

    auto collider_cell_locator = vtkSmartPointer<vtkModifiedBSPTree>::New();
    collider_cell_locator->SetDataSet(collider_polydata);
    collider_cell_locator->Update();

    auto tmp_points_mesh = vtkSmartPointer<vtkPoints>::New();
    auto tmp_points_collider = vtkSmartPointer<vtkPoints>::New();
    auto tmp_cells_collider = vtkSmartPointer<vtkIdList>::New();

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
                contacts.push_back({ .pid = i,
                                     .dx = static_cast<float>(-mesh_normal[0] * distance),
                                     .dy = static_cast<float>(-mesh_normal[1] * distance),
                                     .dz = static_cast<float>(-mesh_normal[2] * distance)});
                //printf("point %d - OK, moving by distance %g\n", i, distance);
            }
        }
    }

    return contacts;
}

void VtkPokeEffect::handle_reaction_simple(vtkPolyData *mesh_polydata, double falloff_radius, double falloff_exponent,
                                           const std::vector<Contact> &contacts, bool debug, double max_distance) {
    auto mesh_normals = mesh_polydata->GetPointData()->GetArray("Normals");

    auto new_mesh_points = vtkSmartPointer<vtkPoints>::New();
    new_mesh_points->DeepCopy(mesh_polydata->GetPoints()); // copy this so iterations do not affect each other

    auto mesh_falloff_locator = vtkSmartPointer<vtkOctreePointLocator>::New();
    mesh_falloff_locator->SetDataSet(mesh_polydata);
    mesh_falloff_locator->Update();

    auto displacement_distance_array = vtkSmartPointer<vtkFloatArray>::New();
    displacement_distance_array->SetNumberOfComponents(1);
    displacement_distance_array->SetNumberOfTuples(mesh_polydata->GetNumberOfPoints());
    displacement_distance_array->FillValue(0.0);

    double min_coord[3] = { DBL_MAX, DBL_MAX, DBL_MAX };
    double max_coord[3] = { DBL_MIN, DBL_MIN, DBL_MIN };

    for (auto c : contacts) {
        float d = std::sqrt(c.dx*c.dx + c.dy*c.dy + c.dz*c.dz);
        displacement_distance_array->SetValue(c.pid, d);

        double p[3];
        mesh_polydata->GetPoint(c.pid, p);
        for (int i = 0; i < 3; i++) {
            min_coord[i] = std::min(min_coord[i], p[i]);
            max_coord[i] = std::max(max_coord[i], p[i]);
        }
    }

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
        double new_distance;

        if (debug) {
            new_distance = old_distance;
        } else {
            new_distance = std::min(max_distance, std::max(old_distance, weighted_average_displacement_distance));
        }

        double new_p[3] = { p0[0] - mesh_normal[0]*new_distance,
                            p0[1] - mesh_normal[1]*new_distance,
                            p0[2] - mesh_normal[2]*new_distance };
        new_mesh_points->SetPoint(i, new_p);
    }

    // final step - overwrite original coordinates with new ones
    mesh_polydata->GetPoints()->ShallowCopy(new_mesh_points);
}
