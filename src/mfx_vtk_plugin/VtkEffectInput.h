#pragma once

#include <vtkPolyData.h>
#include <string>
#include "VtkEffectInputDef.h"

struct VtkEffectInput {
    vtkSmartPointer<vtkPolyData> data;
    double transform[16];
    VtkEffectInputDef *definition;
};
