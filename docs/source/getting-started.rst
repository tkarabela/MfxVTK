Get started with MfxVTK
=======================

To use MfxVTK, you will need the plugin itself as well as a 3D host application
that supports the Open Mesh Effect standard.

Get Open Mesh Effect of Blender
-------------------------------

First, you will need the 3D host. Get the Open Mesh Effect branch of Blender
by `@eliemitchel <https://twitter.com/exppad/>`_ on GitHub: https://github.com/eliemichel/OpenMeshEffectForBlender

Binaries for Windows are available, or you can build it yourself. Just follow official
Blender build guide.

Get the MfxVTK plugin
---------------------

Download the binaries
*********************

This is the simplest way to get started. Just download the ``.ofx`` file
with the MfxVTK Open Mesh Effect plugin from GitHub: https://github.com/tkarabela/MfxVTK/releases

MfxVTK binaries are available for Windows and Linux.

Build MfxVTK yourself
*********************

.. tip::

    Feel free to skip this section if you're happy to use the GitHub binaries.

You will need Git, CMake and a C++ compiler toolchain (like Microsoft Visual Studio or GCC/Make).

First, build VTK (version >= 9.0), see `instructions on VTK wiki <https://vtk.org/Wiki/VTK/Configure_and_Build>`_.

- Use ``BUILD_SHARED_LIBS=OFF`` to make your MfxVTK binary not depend on VTK shared libraries.
- Turn off ``VTK_GROUP_ENABLE_Rendering``, ``VTK_GROUP_ENABLE_Web``, ``VTK_GROUP_ENABLE_MPI``, ``VTK_GROUP_ENABLE_QT`` to make compilation faster.
- (Optional) Use ``VTK_SMP_IMPLEMENTATION_TYPE`` other than ``Sequential`` for better performance of some effects.

Then clone MfxVTK from GitHub (note that it has the OpenMeshEffect repository as a submodule):

.. code-block:: bash

    git clone --recurse-submodules https://github.com/tkarabela/MfxVTK

Build MfxVTK using CMake:

.. code-block:: bash

    cd MfxVTK
    mkdir cmake-build
    cd cmake-build
    cmake .. -DVTK_DIR:PATH=/path/to/VTK-9.0.1-build  # change this accordingly
    cmake --build . --config Debug

The plugin is now in your build directory: ``src/mfx_vtk_plugin/libmfx_vtk_plugin.ofx``.

Use the Open Mesh Effect Modifier
---------------------------------

- Select the default cube
- Go to the Modifier properties tab
- Click *Add modifier...* > *Open Mesh Effect*
- Set *Plugin path* to the MfxVTK ``.ofx`` file you downloaded (or built)
- From the *Select an effect* drop-down menu, pick *Sample points (volume)*
- Et voilà, you're using MfxVTK!

.. tip::
    Try editing the cube in Edit mode. The points will be updated instantly
    to fill the new shape!

.. figure:: /_static/get-started-default-cube.png

    The default cube with *Sample points (volume)* effect, 2000 points, *Distribute uniformly* is off.
    To show the bounding box, turn on *Bounds* in *Viewport display* on the *Object properties* tab.
