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
    const char *PARAM_OFFSET = "Offset";
    const char *PARAM_DEBUG = "Debug";

public:
    struct Contact { int pid; float dx, dy, dz; };

    const char* GetName() override;
    OfxStatus vtkDescribe(OfxParamSetHandle parameters) override;
    OfxStatus vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) override;
    static OfxStatus vtkCook_inner(vtkPolyData *input_polydata, vtkPolyData *output_polydata,
                                   double max_distance, double falloff_radius, double falloff_exponent,
                                   double offset, bool debug);

    static std::vector<Contact> evaluate_collision(vtkPolyData *mesh_polydata, vtkPolyData *collider_polydata,
                                                   double max_distance, double offset);

    static void handle_reaction_simple(vtkPolyData *mesh_polydata, double falloff_radius, double falloff_exponent,
                                       const std::vector<Contact> &contacts, bool debug, double max_distance);
};
