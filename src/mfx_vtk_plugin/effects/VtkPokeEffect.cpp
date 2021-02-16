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
#include "VtkDistanceAlongSurfaceEffect.h"
#include <chrono>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkIdFilter.h>
#include <vtkFloatArray.h>
#include <vtkPolyDataNormals.h>
#include <vtkThreshold.h>
#include <vtkGeometryFilter.h>
#include <vtkModifiedBSPTree.h>
#include <vtkOctreePointLocator.h>
#include <vtkStaticCellLinks.h>
#include <vtkWarpVector.h>
#include <vtkXMLPolyDataWriter.h>

const char *VtkPokeEffect::GetName() {
    return "Poke";
}

OfxStatus VtkPokeEffect::vtkDescribe(OfxParamSetHandle parameters, VtkEffectInputDef &input_mesh, VtkEffectInputDef &output_mesh) {
    input_mesh.RequestVertexAttribute("color0", 3, MfxAttributeType::UByte, MfxAttributeSemantic::Color, true);

    // TODO support collision with other mesh
    // AddInput("Collider");

    // TODO support specifying attribute for using main mesh as collider

    AddParam(PARAM_MAX_DISTANCE, 0.0).Range(0.0, 1e6).Label("Maximum distance");
    AddParam(PARAM_FALLOFF_RADIUS, 0.01).Range(0.0, 1e6).Label("Falloff radius");
    AddParam(PARAM_FALLOFF_EXPONENT, 2.0).Range(1e-2, 1e2).Label("Falloff exponent");
    AddParam(PARAM_COLLISION_SMOOTHING_RATIO, 0.001).Range(0.0, 1.0).Label("Collision smoothing ratio");
    AddParam(PARAM_COLLIDER_NORMAL_FACTOR, 0.0).Range(0.0, 1.0).Label("Collider normal factor");
    AddParam(PARAM_NUMBER_OF_ITERATIONS, 20).Range(1, 10000).Label("Number of iterations");
    AddParam(PARAM_OFFSET, 0.0).Range(0.0, 1e6).Label("Offset");
    AddParam(PARAM_DEBUG, false).Label("Debug");
    return kOfxStatOK;
}

OfxStatus VtkPokeEffect::vtkCook(VtkEffectInput &main_input, VtkEffectInput &main_output, std::vector<VtkEffectInput> &extra_inputs) {
    auto max_distance = GetParam<double>(PARAM_MAX_DISTANCE).GetValue();
    auto falloff_radius = GetParam<double>(PARAM_FALLOFF_RADIUS).GetValue();
    auto falloff_exponent = GetParam<double>(PARAM_FALLOFF_EXPONENT).GetValue();
    auto collision_smoothing_ratio = GetParam<double>(PARAM_COLLISION_SMOOTHING_RATIO).GetValue();
    auto collider_normal_factor = GetParam<double>(PARAM_COLLIDER_NORMAL_FACTOR).GetValue();
    auto number_of_iterations = GetParam<int>(PARAM_NUMBER_OF_ITERATIONS).GetValue();
    auto offset = GetParam<double>(PARAM_OFFSET).GetValue();
    auto debug = GetParam<bool>(PARAM_DEBUG).GetValue();

    if (!main_input.data->GetPointData()->HasArray("color0")) {
        printf("VtkPokeEffect - error, input must have vertex color attribute called 'color0'\n");
        return kOfxStatFailed;
    }

    // XXX NOTE TO USER - paint mesh white and collider black...

    vtkCook_inner(main_input.data, main_output.data, max_distance, falloff_radius, falloff_exponent,
                  collision_smoothing_ratio, offset, number_of_iterations, debug, collider_normal_factor);

    return kOfxStatOK;
}

OfxStatus VtkPokeEffect::vtkCook_inner(vtkPolyData *input_polydata, vtkPolyData *output_polydata, double max_distance,
                                       double falloff_radius,
                                       double falloff_exponent, double collision_smoothing_ratio, double offset,
                                       int number_of_iterations,
                                       bool debug, double collider_normal_factor) {
    auto t0 = std::chrono::system_clock::now();

    // remember point IDs so that we can map deformed mesh onto input_polydata
    // TODO even better would be to keep input vtkPoints as much as possible
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

    auto contacts = evaluate_collision(mesh_polydata, collider_polydata, max_distance, offset, debug, collider_normal_factor);
    auto t2 = std::chrono::system_clock::now();
    handle_reaction_laplacian(mesh_polydata, contacts, falloff_radius, falloff_exponent, number_of_iterations, collision_smoothing_ratio);
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
                                  double offset, bool debug, double collider_normal_factor) {
    std::vector<Contact> contacts;
    int n = mesh_polydata->GetNumberOfPoints();
    double mesh_diagonal_length = mesh_polydata->GetLength();
    
    vtkFloatArray *this_to_collider_arr = nullptr, *this_to_mesh_arr = nullptr, *mesh_to_collider_arr = nullptr, *contact_arr = nullptr;
    
    if (debug) {
        this_to_collider_arr = vtkFloatArray::New();
        this_to_collider_arr->SetNumberOfComponents(1);
        this_to_collider_arr->SetNumberOfTuples(n);
        this_to_collider_arr->FillValue(vtkMath::Nan());
        this_to_collider_arr->SetName("this_to_collider");
        mesh_polydata->GetPointData()->AddArray(this_to_collider_arr);

        this_to_mesh_arr = vtkFloatArray::New();
        this_to_mesh_arr->SetNumberOfComponents(1);
        this_to_mesh_arr->SetNumberOfTuples(n);
        this_to_mesh_arr->FillValue(vtkMath::Nan());
        this_to_mesh_arr->SetName("this_to_mesh");
        mesh_polydata->GetPointData()->AddArray(this_to_mesh_arr);

        mesh_to_collider_arr = vtkFloatArray::New();
        mesh_to_collider_arr->SetNumberOfComponents(1);
        mesh_to_collider_arr->SetNumberOfTuples(n);
        mesh_to_collider_arr->FillValue(vtkMath::Nan());
        mesh_to_collider_arr->SetName("mesh_to_collider");
        mesh_polydata->GetPointData()->AddArray(mesh_to_collider_arr);

        contact_arr = vtkFloatArray::New();
        contact_arr->SetNumberOfComponents(3);
        contact_arr->SetNumberOfTuples(n);
        contact_arr->FillValue(0);
        contact_arr->SetName("contact");
        mesh_polydata->GetPointData()->AddArray(contact_arr);
    }

    if (offset > 0) {
        auto warp_filter = vtkSmartPointer<vtkWarpVector>::New();
        warp_filter->SetInputData(collider_polydata);
        warp_filter->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "Normals");
        warp_filter->SetScaleFactor(offset);
        warp_filter->Update();
        collider_polydata->SetPoints(warp_filter->GetOutput()->GetPoints());
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
    for (int i = 0; i < n; i++) {
        double p[3];
        double mesh_point[3], collider_point[3];
        double mesh_normal[3], collider_cell_normal[3];
        mesh_polydata->GetPoint(i, p);

        // early exit, mesh point too far away from collider
        if (p[0] < min_coord_estimate[0] || p[0] > max_coord_estimate[0] ||
            p[1] < min_coord_estimate[1] || p[1] > max_coord_estimate[1] ||
            p[2] < min_coord_estimate[2] || p[2] > max_coord_estimate[2]) {
            continue;
        }

        mesh_normals->GetTuple(i, mesh_normal);

        // Raycasting setup, where each point will fire a ray "behind" it in normal direction.
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

        if (debug) {
            this_to_collider_arr->SetValue(i, std::sqrt(this_to_collider_sq));
            this_to_mesh_arr->SetValue(i, std::sqrt(this_to_mesh_sq));
            mesh_to_collider_arr->SetValue(i, std::sqrt(mesh_to_collider_sq));
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
                double v[3];
                for (int j = 0; j < 3; j++) {
                    v[j] = -mesh_normal[j]*(1.0-collider_normal_factor) + collider_cell_normal[j]*collider_normal_factor;
                }

                Contact contact = { .pid = i,
                                    .dx = static_cast<float>(v[0] * distance),
                                    .dy = static_cast<float>(v[1] * distance),
                                    .dz = static_cast<float>(v[2] * distance)};
                contacts.push_back(contact);

                if (debug) {
                    contact_arr->SetTuple3(i, contact.dx, contact.dy, contact.dz);
                }
                //printf("point %d - OK, moving by distance %g\n", i, distance);
            }
        }
    }

    if (debug) {
        {
            auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
            writer->SetInputData(mesh_polydata);
            printf("VtkPokeEffect - dumping mesh to vtkpokeeffect-mesh.vtp\n");
            writer->SetFileName("vtkpokeeffect.vtp");
            writer->Write();
        }
        {
            auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
            writer->SetInputData(collider_polydata);
            printf("VtkPokeEffect - dumping collider to vtkpokeeffect-collider.vtp\n");
            writer->SetFileName("vtkpokeeffect-collider.vtp");
            writer->Write();
        }
    }

    return contacts;
}

void VtkPokeEffect::handle_reaction_laplacian(vtkPolyData *mesh_polydata, const std::vector<Contact> &contacts,
                                              double falloff_radius, double falloff_exponent,
                                              int number_of_iterations, double collision_smoothing_ratio) {
    auto mesh_normals = mesh_polydata->GetPointData()->GetArray("Normals");

    // apply initial deformation
    auto new_mesh_points = vtkSmartPointer<vtkPoints>::New();
    new_mesh_points->DeepCopy(mesh_polydata->GetPoints());

    std::vector<int> collision_points;
    for (auto c : contacts) {
        collision_points.push_back(c.pid);

        double p[3];
        new_mesh_points->GetPoint(c.pid, p);
        p[0] += c.dx;
        p[1] += c.dy;
        p[2] += c.dz;
        new_mesh_points->SetPoint(c.pid, p);
    }

    auto cell_links = vtkSmartPointer<vtkStaticCellLinks>::New(); // TODO use templated class to get int32 here
    cell_links->BuildLinks(mesh_polydata);

    auto manifold_distance_arr = VtkDistanceAlongSurfaceEffect::compute_distance(mesh_polydata, collision_points.size(),
                                                                                 collision_points.data(), cell_links,
                                                                                 (float) falloff_radius);

    auto tmp_id_list = vtkSmartPointer<vtkIdList>::New(); // TODO get rid of this, prevents parallelization
    auto push_point_neighbors = [&mesh_polydata, &cell_links, &tmp_id_list](int u, std::vector<int> &connectivity) -> int {
        tmp_id_list->Reset();
        int num_cells = cell_links->GetNumberOfCells(u);
        auto neighbor_cells = cell_links->GetCells(u);
        for (int i = 0; i < num_cells; i++) {
            auto cell = mesh_polydata->GetCell(neighbor_cells[i]);
            int num_points = cell->GetNumberOfPoints();

            for (int j = 0; j < num_points; j++) {
                int v = cell->GetPointId(j);
                tmp_id_list->InsertUniqueId(v);
            }
        }

        for (int i = 0; i < tmp_id_list->GetNumberOfIds(); i++) {
            connectivity.push_back(tmp_id_list->GetId(i));
        }
        return tmp_id_list->GetNumberOfIds();
    };

    // pick points to be affected by the laplacian
    std::vector<int> smoothed_points;
    std::vector<int> smoothed_points_offsets;
    std::vector<int> smoothed_points_connectivity;
    int previous_offset = 0;
    for (int i = 0; i < manifold_distance_arr->GetNumberOfValues(); i++) {
        float d = manifold_distance_arr->GetValue(i);

        // smooth points in falloff_radius
        // if collision_smoothing_ratio > 0, we want to smooth collision points as well
        if (d < falloff_radius && (d > 0 || collision_smoothing_ratio > 0)) {
            smoothed_points.push_back(i);
            smoothed_points_offsets.push_back(previous_offset);
            int neighbor_count = push_point_neighbors(i, smoothed_points_connectivity);
            previous_offset += neighbor_count;
        }
    }
    smoothed_points_offsets.push_back(previous_offset);
    printf("VtkPokeEffect - picked %d points for smoothing\n", smoothed_points.size());

    // iterate laplacian
    for (int i = 0; i < number_of_iterations; i++) {
        // TODO parallelize
        // XXX make iterations independent!!!
        for (int j = 0; j < smoothed_points.size(); j++) {
            int pid0 = smoothed_points[j];
            double p0[3];
            new_mesh_points->GetPoint(pid0, p0);
            double barycenter[3] = {0, 0, 0};
            int num_neighbors = smoothed_points_offsets[j+1] - smoothed_points_offsets[j];

            for (int k = smoothed_points_offsets[j]; k < smoothed_points_offsets[j+1]; k++) {
                int pid = smoothed_points_connectivity[k];
                double p[3];
                new_mesh_points->GetPoint(pid, p);
                barycenter[0] += p[0];
                barycenter[1] += p[1];
                barycenter[2] += p[2];
            }
            barycenter[0] /= num_neighbors;
            barycenter[1] /= num_neighbors;
            barycenter[2] /= num_neighbors;

            float distance_to_collision = manifold_distance_arr->GetValue(pid0);
            double alpha = 0.0;  // blending weight

            double disp[3] = { barycenter[0] - p0[0],
                               barycenter[1] - p0[1],
                               barycenter[2] - p0[2] };

            if (distance_to_collision > 0) {
                alpha = std::pow(std::max(0.0, 1.0 - (distance_to_collision / falloff_radius)), falloff_exponent);
            } else {
                // collision point
                double normal[3];
                mesh_normals->GetTuple(pid0, normal);
                if (vec3_dot(disp, normal) > 0) {
                    // pushing outside - back towards collider, we want to dampen this
                    alpha = collision_smoothing_ratio;
                } else {
                    // pushing inside mesh, use full laplacian
                    alpha = 1.0;
                }
            }

            p0[0] += alpha * disp[0];
            p0[1] += alpha * disp[1];
            p0[2] += alpha * disp[2];

            new_mesh_points->SetPoint(pid0, p0);
        }
    }

    // final step - overwrite original coordinates with new ones
    mesh_polydata->GetPoints()->ShallowCopy(new_mesh_points);
}
