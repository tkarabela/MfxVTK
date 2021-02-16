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

#pragma once

#include "VtkEffect.h"

class VtkPokeEffect : public VtkEffect {
private:
    const char *PARAM_MAX_DISTANCE = "MaxDistance";
    const char *PARAM_FALLOFF_RADIUS = "FalloffRadius";
    const char *PARAM_FALLOFF_EXPONENT = "FalloffExponent";
    const char *PARAM_COLLISION_SMOOTHING_RATIO = "CollisionSmoothingRatio";
    const char *PARAM_COLLIDER_NORMAL_FACTOR = "ColliderNormalFactor";
    const char *PARAM_NUMBER_OF_ITERATIONS = "NumberOfIterations";
    const char *PARAM_OFFSET = "Offset";
    const char *PARAM_DEBUG = "Debug";

public:
    struct Contact { int pid; float dx, dy, dz; };

    const char* GetName() override;
    OfxStatus vtkDescribe(OfxParamSetHandle parameters, VtkEffectInputDef &input_mesh, VtkEffectInputDef &output_mesh) override;
    OfxStatus vtkCook(VtkEffectInput &main_input, VtkEffectInput &main_output, std::vector<VtkEffectInput> &extra_inputs) override;
    static OfxStatus
    vtkCook_inner(vtkPolyData *input_polydata, vtkPolyData *output_polydata, double max_distance, double falloff_radius,
                  double falloff_exponent, double collision_smoothing_ratio, double offset, int number_of_iterations,
                  bool debug, double collider_normal_factor);

    /* Find out which mesh points need to be moved to clear the collision.
     * */
    static std::vector<Contact>
    evaluate_collision(vtkPolyData *mesh_polydata, vtkPolyData *collider_polydata, double max_distance, double offset,
                       bool debug, double collider_normal_factor);

    /* Clear collision by depressing points along their normals;
     * create falloff by laplacian smoothing weighted by manifold distance from the collision.
     * */
    static void handle_reaction_laplacian(vtkPolyData *mesh_polydata, const std::vector<Contact> &contacts,
                                          double falloff_radius, double falloff_exponent, int number_of_iterations,
                                          double collision_smoothing_ratio);
};
