Sample Points (Surface)
***********************

Generate points along edges or faces at regular distance.

.. note::
    Random sampling should be coming in next VTK release (`reference <https://github.com/Kitware/VTK/commit/c246e3dd3e28b513df521b2ccfe2a34bb83a6d2a>`_).

:Input: edge wireframe, polygonal mesh
:Output: point cloud
:VTK classes: ``vtkPolyDataPointSampler``, ``vtkTriangleFilter``

Options
#######

.. figure:: /_static/properties-sample-points-surface.png
    :align: right
    :width: 300px

Distance
    Target spacing of the sampled points.

Sample edges
    If this option is turned on, points will be sampled along edges.

Sample faces
    If this option is turned on, points will be sampled on faces.

Example
#######

.. figure:: /_static/example-sample-points-surface.png

    Left: The default cube.
    Middle: Cube with *Sample points (surface)* and *Sample edges*.
    Right: Cube with *Sample points (surface)* and *Sample faces*.
