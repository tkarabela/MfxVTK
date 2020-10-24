Decimate
********

This effect will reduce the number of faces of a mesh while trying to preserve its shape.

.. note::
    This will triangulate the mesh first.

:Input: polygonal mesh
:Output: polygonal mesh
:VTK classes: ``vtkQuadricDecimation``, ``vtkTriangleFilter``

Options
#######

.. figure:: /_static/properties-decimate.png
    :align: right
    :width: 300px


Target ratio
    Ratio of input : output triangles. For example, *Target ratio* = 0.2 means 80% reduction
    of triangle count, with the decimated mesh only having 20% triangles compared to the original.

Preserve volume
    When this option is on, mesh shape will be more accurately preserved, at cost of performance.

Example
#######

.. figure:: /_static/example-decimate.png

    Left: The Stanford bunny.
    Middle: The bunny with *Decimate* effect (*Target ratio* = 0.1).
    Right: The bunny with *Decimate* effect (*Target ratio* = 0.01).
