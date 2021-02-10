#include <vtkMatrix4x4.h>
#include "VtkEffectInput.h"


void VtkEffectInput::get_relative_transform(const VtkEffectInput &other, double *output_transorm) const {
    double m1_inv[16];
    vtkMatrix4x4::Invert(this->transform, m1_inv);
    vtkMatrix4x4::Multiply4x4(m1_inv, other.transform, output_transorm);
}
