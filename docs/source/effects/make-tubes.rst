Make Tubes
**********

This effect will create *n*-sided cylinders along given edges.

.. tip::
   This is especially helpful if your input doesn't have polygons.

   Blender's *Wireframe* modifier creates nicer connections where edges meet,
   but it requires polygonal mesh, making it less composable with other effects.

:Input: edge wireframe, polygonal mesh
:Output: polygonal mesh
:VTK classes: ``vtkTubeFilter``, ``vtkExtractEdges``

Options
#######

.. figure:: /_static/properties-make-tubes.png
    :align: right
    :width: 300px

Radius
    Radius of the cylinder created along each edge.

Number of sides
    Number of sides for the cylinders (more sides create nicer looking mesh,
    but are more resource intensive).

Capping
    Whether the ends should be closed.

Example
#######

.. figure:: /_static/example-make-tubes.png

    Left: Edges of an icosphere. Right: Icosphere edges with the "Make tubes" effect.
