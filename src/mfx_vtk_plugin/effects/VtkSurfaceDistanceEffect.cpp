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

#include "VtkSurfaceDistanceEffect.h"
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vector>
#include <vtkStaticCellLinks.h>
#include <queue>

const char *VtkSurfaceDistanceEffect::GetName() {
    return "Surface distance";
}

OfxStatus VtkSurfaceDistanceEffect::vtkDescribe(OfxParamSetHandle parameters) {
    // TODO add input color attribute when it lands in OFX
    // TODO declare this is a deformer

    AddParam(PARAM_NORMALIZE_DISTANCE, true).Label("Normalize distance");
    return kOfxStatOK;
}

OfxStatus VtkSurfaceDistanceEffect::vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) {
    auto input_color_arr = input_polydata->GetPointData()->GetArray("color0");
    auto normalize_distance = GetParam<bool>(PARAM_NORMALIZE_DISTANCE).GetValue();

    if (!input_color_arr) {
        printf("VtkSurfaceDistanceEffect - missing input 'color0' attribute!\n");
        return kOfxStatFailed;
    }

    int n = input_polydata->GetNumberOfPoints();
    std::vector<int> source_points;

    for (int i = 0; i < n; i++) {
        if (input_color_arr->GetComponent(i, 0) > 0) {
            // has nonzero in Red channel
            source_points.push_back(i);
        }
    }

    if (source_points.empty()) {
        printf("VtkSurfaceDistanceEffect - warning, no source points!\n");
    } else {
        printf("VtkSurfaceDistanceEffect - I have %d source points\n", source_points.size());
    }

    auto cell_links = vtkSmartPointer<vtkStaticCellLinks>::New(); // TODO use templated class to get int32 here
    cell_links->BuildLinks(input_polydata);

    auto distance_arr = compute_distance(input_polydata, source_points.size(), source_points.data(), cell_links);

    auto output_uv_arr = vtkFloatArray::New();
    output_uv_arr->SetNumberOfComponents(2);
    output_uv_arr->SetNumberOfTuples(n);

    if (normalize_distance) {
        float range[2];
        distance_arr->GetFiniteValueRange(range);
        float max_distance = range[1];
        printf("VtkSurfaceDistanceEffect - normalizing by max value %g\n", max_distance);
        for (int i = 0; i < n; i++) {
            float value = distance_arr->GetValue(i);
            if (!vtkMath::IsFinite(value)) {
                value = 0.0f;
            }
            float norm_distance = value / max_distance;
            output_uv_arr->SetTypedComponent(i, 0, norm_distance);
        }
    } else {
        // TODO handle inf here
        output_uv_arr->CopyComponent(0, distance_arr, 0);
    }

    output_uv_arr->FillTypedComponent(1, 0.0);
    output_uv_arr->SetName("uv0"); // TODO handle output better w.r.t. existing UV arrays

    output_polydata->ShallowCopy(input_polydata);
    output_polydata->GetPointData()->AddArray(output_uv_arr);

    return kOfxStatOK;
}

vtkFloatArray *VtkSurfaceDistanceEffect::compute_distance(vtkPolyData *mesh, int num_source_points,
                                                          const int *source_points,
                                                          vtkStaticCellLinks *cell_links, float max_distance) {
    int n = mesh->GetNumberOfPoints();

    auto manifold_distance_arr = vtkFloatArray::New();
    manifold_distance_arr->SetNumberOfComponents(1);
    manifold_distance_arr->SetNumberOfTuples(n);
    manifold_distance_arr->SetName("ManifoldDistance");
    manifold_distance_arr->FillValue(vtkMath::Inf());

    enum class Status : unsigned char { UNVISITED, OPEN, CLOSED, SOURCE };
    std::vector<Status> status_arr(n, Status::UNVISITED);

    // add source points to queue
    typedef std::pair<float, int> PointDistance;
    auto cmp = [](PointDistance left, PointDistance right) { return left.first > right.first; };
    std::priority_queue<PointDistance, std::vector<PointDistance>, decltype(cmp)> queue(cmp);

    for (int i = 0; i < num_source_points; i++) {
        int p = source_points[i];
        status_arr[p] = Status::SOURCE; // for debugging only
        manifold_distance_arr->SetValue(p, 0.0);
        queue.push(std::make_pair(0.0, p));
    }

    auto get_distance = [&mesh](int u, int v) -> float {
        double x[3], y[3];
        mesh->GetPoints()->GetPoint(u, x);
        mesh->GetPoints()->GetPoint(v, y);
        double dist_sq = vtkMath::Distance2BetweenPoints(x, y);
        return (float)std::sqrt(dist_sq);
    };

    // begin dijkstra
    int it = 0;
    while (!queue.empty()) {
        float u_dist; int u;
        std::tie(u_dist, u) = queue.top();
        queue.pop();
        it++;

        if (u_dist > max_distance) {
            // early exit, we've computed all closest paths up to max_distance
            break;
        }

        Status u_status = status_arr[u];
        float u_distance = manifold_distance_arr->GetValue(u);
        //printf("iteration %d - point %d, dist %g, status %d\n", it, u, u_dist, (int)u_status);

        if (u_status == Status::CLOSED) {
            continue;
        };

        int num_cells = cell_links->GetNumberOfCells(u);
        auto neighbor_cells = cell_links->GetCells(u);
        for (int i = 0; i < num_cells; i++) {
            auto cell = mesh->GetCell(neighbor_cells[i]);
            int num_points = cell->GetNumberOfPoints();

            for (int j = 0; j < num_points; j++) {
                int v = cell->GetPointId(j);
                Status v_status = status_arr[v];

                float uv_distance = get_distance(u, v);
                float old_v_distance = manifold_distance_arr->GetValue(v);
                float new_v_distance = u_distance + uv_distance;

                //printf("  - point %d, old dist %g, new dist %g, status %d\n", v, old_v_distance, new_v_distance, (int)v_status);

                if (new_v_distance < old_v_distance) {
                    manifold_distance_arr->SetValue(v, new_v_distance);
                    queue.push(std::make_pair(new_v_distance, v));
                    status_arr[v] = Status::OPEN;

                    //printf("  - enqueue point %d\n", v);
                }
            }
            // end neighbor cell
        }
        status_arr[u] = Status::CLOSED;
        // end point
    }
    // end dijkstra

    return manifold_distance_arr;
}
