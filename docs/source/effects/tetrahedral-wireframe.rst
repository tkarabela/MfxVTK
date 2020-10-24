Tetrahedral Wireframe
*********************

This effect creates a wireframe from a point cloud. It is based on Delaunay triangulation
in 3D, which creates a convex hull filled with tetrahedra. To reveal concave features, use
the *Maximum length* parameter, which will eat long edges from the convex hull.

.. note::
    In the future, there may be a more proper option to handle this (alpha-shapes,
    concave hull).

:Input: point cloud
:Output: edge wireframe
:VTK classes: ``vtkDelaunay3D``, ``vtkExtractEdges``, ``vtkCleanPolyData``

Options
#######

.. figure:: /_static/properties-tetrahedral-wireframe.png
    :align: right
    :width: 300px


Maximum edge length
    This specifies maximum length for edges. Setting a lower value will create
    a more defined shape. Too small values may create blocky or disconnected mesh.


Example
#######

.. figure:: /_static/example-tetrahedral-wireframe.png

    Left: Stanford bunny without modifiers. Right: Stanford bunny with "Sample points (volume)"
    at 1000 points, then "Tetrahedral wireframe" effect.