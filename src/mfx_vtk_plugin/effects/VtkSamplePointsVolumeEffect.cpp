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

#include <vtkFloatArray.h>
#include <vtkQuadricDecimation.h>
#include <vtkMinimalStandardRandomSequence.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkTriangleFilter.h>
#include <vtkPointData.h>

#include "VtkSamplePointsVolumeEffect.h"
#include "mfx_vtk_utils.h"

const char *VtkSamplePointsVolumeEffect::GetName() {
    return "Sample points (volume)";
}

OfxStatus VtkSamplePointsVolumeEffect::vtkDescribe(OfxParamSetHandle parameters, VtkEffectInputDef &input_mesh,
                                                   VtkEffectInputDef &output_mesh) {
    AddParam(PARAM_NUMBER_OF_POINTS, 200).Range(1, 1e6).Label("Number of points");
    AddParam(PARAM_DISTRIBUTE_UNIFORMLY, true).Label("Distribute points uniformly");
    AddParam(PARAM_AUTO_SIMPLIFY, true).Label("Auto simplify input mesh");
    // TODO more controls
    return kOfxStatOK;
}

OfxStatus VtkSamplePointsVolumeEffect::vtkCook(VtkEffectInput &main_input, VtkEffectInput &main_output, std::vector<VtkEffectInput> &extra_inputs) {
    auto number_of_points = GetParam<int>(PARAM_NUMBER_OF_POINTS).GetValue();
    auto distribute_uniformly = GetParam<bool>(PARAM_DISTRIBUTE_UNIFORMLY).GetValue();
    auto auto_simplify = GetParam<bool>(PARAM_AUTO_SIMPLIFY).GetValue();
    return vtkCook_inner(main_input.data, main_output.data, number_of_points, distribute_uniformly, auto_simplify);
}

OfxStatus VtkSamplePointsVolumeEffect::vtkCook_inner(vtkPolyData *input_polydata, vtkPolyData *output_polydata,
                                                     int number_of_points, bool distribute_uniformly,
                                                     bool auto_simplify, bool _assume_input_polydata_triangles) {
    double bounds[6];
    input_polydata->GetBounds(bounds);

    auto poly_data_distance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();

    if (auto_simplify && input_polydata->GetNumberOfPolys() > 100) {
        vtkSmartPointer<vtkPolyData> input_triangle_mesh = input_polydata;
        if (!_assume_input_polydata_triangles) {
            auto triangle_filter = vtkSmartPointer<vtkTriangleFilter>::New();
            triangle_filter->SetInputData(input_polydata);
            triangle_filter->Update();
            input_triangle_mesh = triangle_filter->GetOutput();
        }

        auto input_polycount = input_triangle_mesh->GetNumberOfPolys();
        int target_polycount = 1000 + static_cast<int>(std::sqrt(input_polycount));
        double target_reduction = static_cast<double>(target_polycount) / input_polycount;

        if (target_reduction < 1) {
            auto quadratic_decimation_filter = vtkSmartPointer<vtkQuadricDecimation>::New();
            quadratic_decimation_filter->SetInputData(input_triangle_mesh);
            quadratic_decimation_filter->SetTargetReduction(target_reduction);
            quadratic_decimation_filter->Update();
            poly_data_distance->SetInput(quadratic_decimation_filter->GetOutput());
        } else {
            poly_data_distance->SetInput(input_polydata);
        }
    } else {
        poly_data_distance->SetInput(input_polydata);
    }

    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(number_of_points);

    output_polydata->SetPoints(points);

    auto distance_arr = vtkSmartPointer<vtkFloatArray>::New();
    distance_arr->SetName("distance");
    distance_arr->SetNumberOfComponents(1);
    distance_arr->SetNumberOfTuples(number_of_points);

    output_polydata->GetPointData()->AddArray(distance_arr);

    auto random_generator_vtk = vtkSmartPointer<vtkMinimalStandardRandomSequence>::New();
    auto random_generator_custom = AdditiveRecurrence<3>();
    auto get_random_uniform = [distribute_uniformly, &random_generator_vtk, &random_generator_custom](int i, double low, double high) -> double {
        if (distribute_uniformly) {
            random_generator_custom.Next();
            return random_generator_custom.GetRangeValue(i, low, high);
        } else {
            random_generator_vtk->Next();
            return random_generator_vtk->GetRangeValue(low, high);
        };
    };

    int i, iteration_count;
    for (i = 0, iteration_count = 0;
         i < number_of_points && iteration_count < 10*number_of_points;
         iteration_count++)
    {
        double x = get_random_uniform(0, bounds[0], bounds[1]),
               y = get_random_uniform(1, bounds[2], bounds[3]),
               z = get_random_uniform(2, bounds[4], bounds[5]);

        double distance = poly_data_distance->EvaluateFunction(x, y, z);

        if (distance < 0) {
            points->SetPoint(i, x, y, z);
            distance_arr->SetValue(i, distance);
            i++;
        }
    }

    if (i < number_of_points) {
        printf("WARNING - gave up after %d iterations, but I only have %d points\n", iteration_count, i);
        points->SetNumberOfPoints(i);
        distance_arr->SetNumberOfTuples(i);
    }

    // TODO write error if input has no faces
    // TODO use output distance array

    return kOfxStatOK;
}
