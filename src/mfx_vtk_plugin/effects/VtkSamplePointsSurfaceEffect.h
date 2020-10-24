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

class VtkSamplePointsSurfaceEffect : public VtkEffect {
private:
    const char *PARAM_DISTANCE = "Distance";
    // const char *PARAM_DISTRIBUTE_UNIFORMLY = "DistributeUniformly"; // TODO enable this when it lands in VTK >9.0
    // const char *PARAM_GENERATE_VERTEX_POINTS = "GenerateVertexPoints";
    const char *PARAM_GENERATE_EDGE_POINTS = "GenerateEdgePoints";
    const char *PARAM_GENERATE_INTERIOR_POINTS = "GenerateInteriorPoints";
    //const char *PARAM_INTERPOLATE_POINT_DATA = "InterpolatePointData";

public:
    const char* GetName() override;
    OfxStatus vtkDescribe(OfxParamSetHandle parameters) override;
    OfxStatus vtkCook(vtkPolyData *input_polydata, vtkPolyData *output_polydata) override;
    static OfxStatus vtkCook_inner(vtkPolyData *input_polydata, vtkPolyData *output_polydata,
                                   double distance, bool generate_vertex_points,
                                   bool generate_edge_points, bool generate_interior_points);
};
