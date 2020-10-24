Reduce Edges
************

This effect removes long edges from the mesh based on a threshold.

.. tip::
    This effect is very handy together with :doc:`Tetrahedral wireframe`,
    to make the resulting mesh more well-defined.

:Input: edge wireframe
:Output: edge wireframe
:VTK classes: ``vtkCleanPolyData``

Options
#######

.. figure:: /_static/properties-reduce-edges.png
    :align: right
    :width: 300px

Maximum length
    This specifies maximum length for edges. Any longer edges will be removed.

Example
#######

.. figure:: /_static/example-reduce-edges.png

    Left: Stanford bunny with "Sample points (volume)" and "Tetrahedral wireframe" effects.
    Right: The same bunny with additional "Reduce edges" effect applied. Notice that the shape is
    much more well-defined, long edges from the convex hull are removed.
