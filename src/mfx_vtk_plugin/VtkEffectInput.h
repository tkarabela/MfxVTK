#pragma once

#include <vtkPolyData.h>
#include <string>
#include "VtkEffectInputDef.h"

struct VtkEffectInput {
    vtkSmartPointer<vtkPolyData> data;
    double transform[16];
    VtkEffectInputDef *definition;

    /*
     * Get transformation matrix so that "vtkTransformFilter(other.data, output_transform)" (pseudocode)
     * is in the same coordinate system as this->data.
     */
    void get_relative_transform(const VtkEffectInput &other, double output_transorm[16]) const;
};
