Comparison to VTK and BVTKNodes
===============================

MfxVTK is not the only way to use The Visualization Toolkit library in Blender.
On this page, you can find comparison to a popular add-on for this purpose, BVTKNodes.

.. |br| raw:: html

   <br />

About VTK
---------

`The Visualization Toolkit <https://vtk.org/>`_ is an open-source data processing and
visualization C++ library used in the sciences, engineering and medicine. Together with related
projects like `Paraview <https://www.paraview.org/>`_ and `ITK <https://itk.org/>`_,
it provides algorithms to handle various types of data, scaling to huge datasets in
High Performance Computing settings, running VR/CAVE sessions, or even interactive
visualizations in the browser via `vtk.js <https://kitware.github.io/vtk-js/>`_.


The Case for BVTKNodes
----------------------

.. figure:: /_static/bvtknodes_motorbike102_cam02.jpg
    :align: right
    :width: 350px

    OpenFOAM motorBike tutorial case visualized |br|
    in Blender with BVTKNodes. Image by `@tkesita <https://github.com/tkeskita>`_.

The `BVTKNodes <https://bvtknodes.readthedocs.io/en/latest/index.html>`_ project
provides a way to create VTK data processing pipelines inside Blender.
This makes it possible to import scientific data into Blender and create photorealistic
visualizations -- a more rendering-oriented alternative to Paraview.

The add-on directly exposes VTK classes using a node editor, which is separate from
rest of Blender. There is a special node to export geometry to Blender mesh for rendering.

In summary, BVTKNodes is specifically aimed at scientific and engineering use (eg. FEM, CFD).


The Case for MfxVTK
-------------------

The MfxVTK library offers new mesh Modifiers through the Open Mesh Effect interface.
It works seamlessly with the rest of Blender, just like native Modifiers. Some of
the effects offer alternatives to stock Modifiers (eg. Decimate, Smooth), while others
bring entirely new functions to the Modifier stack.

Though the Effects are leveraging VTK library behind the scenes, there is no 1:1 mapping
to specific VTK classes. Instead, the Effects are organized with general 3D graphics work in mind
(as opposed to data visualization/processing).

In summary, MfxVTK is aimed at artists and Blender users interested in non-destructive
modelling, without any particular scientific bent.
