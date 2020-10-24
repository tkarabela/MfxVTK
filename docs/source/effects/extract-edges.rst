Extract Edges
*************

This effect can extract edges from mesh based on their properties.

:Input: polygonal mesh
:Output: edge wireframe
:VTK classes: ``vtkFeatureEdges``, ``vtkTriangleFilter``, ``vtkCleanPolydata``

Options
#######

.. figure:: /_static/properties-extract-edges.png
    :align: right
    :width: 300px

Extract feature edges
    Whether "sharp" edges should be extracted.

    .. tip::
        This refers to angle between faces, not explicitly marking edges
        as *Sharp* in Blender.

Feature edge angle
    Minimum angle between faces to consider an edge as a "feature" edge.

Extract boundary edges
    Whether boundary (ie. having just one face) edges should be extracted.

Extract non-manifold edges
    Whether non-manifold (ie. having three or more faces) edges should be extracted.

Extract manifold edges
    Whether manifold (ie. having two faces) edges should be extracted.

Example
#######

.. figure:: /_static/example-extract-edges.png

    Left: Suzanne without modifiers. Right: Suzanne with the "Extract feature edges" turned on.
