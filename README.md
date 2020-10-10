# MfxVTK

*An [Open Mesh Effect][OpenMeshEffect] plugin for [VTK][VTK] filters*

![Smoothing Suzanne the monkey with VTK filters](doc/monkeys.png)

With this plugin, you can process 3D mesh data using filters
from [Visualization Toolkit][VTK] (VTK) library in any software that supports
the [Open Mesh Effect][OpenMeshEffect] standard. For example, they can be used
as modifiers in the [OpenMeshEffectForBlender branch of Blender][OpenMeshEffectForBlender].

🚧 *Currently this is a proof-of-concept; it works, but more features and polish
should come in the future.*

### Features

Effect | VTK filter
------------ | -------------
Laplacian smoothing | `vtkSmoothPolyDataFilter`
Windowed sinc smoothing | `vtkWindowedSincPolyDataFilter`

### How to use it

In the [OpenMeshEffectForBlender branch of Blender][OpenMeshEffectForBlender]:

- select your mesh
- open the **Modifier** tab and select **Open Mesh Effect** modifier
- in modifier properties, set path to MfxVTK file (`mfx_vtk_plugin.ofx`)
- select effect from the drop-down menu

### How to build it yourself

It's a CMake project with dependency on the VTK library (>=9.0). If you're
building VTK yourself, you can point CMake to your build directory, eg.
`-DVTK_DIR:PATH=/path/to/VTK-9.0.1-build`.

🚧 *TODO better explanation*

### How it works

It converts Open Mesh Effect mesh into `vtkPolyData`. General structure
of the plugin as well as handy C++ wrapper over the Open Mesh Effect interface
is borrowed from [MfxVCG].

### License

MfxVTK source code (`src/mfx_vtk_plugin`) is released under the MIT license, see [LICENSE.txt](LICENSE.txt).
Parts of this repository also contain code licensed under the Apache 2 license, please
see copyright notices in individual files.

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

[OpenMeshEffect]: http://openmesheffect.org
[OpenMeshEffectForBlender]: https://github.com/eliemichel/OpenMeshEffectForBlender
[VTK]: https://vtk.org
[MfxVCG]: https://github.com/eliemichel/MfxVCG
