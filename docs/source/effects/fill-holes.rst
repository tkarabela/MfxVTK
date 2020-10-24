Fill Holes
**********

This effect will triangulate holes (loops of boundary edges), which can be handy
for fixing imperfect meshes like 3D scans. This effect only makes sense for an approximately
closed mesh; otherwise the results may be unexpected (eg. if you give it a plane,
it would try to connect the boundaries).

:Input: polygonal mesh (should be approximately closed)
:Output: polygonal mesh
:VTK classes: ``vtkFillHolesFilter``

Options
#######

.. figure:: /_static/properties-fill-holes.png
    :align: right
    :width: 300px

Maximum hole size
    Only fill holes that fit inside a sphere of given radius.

Example
#######

.. figure:: /_static/example-fill-holes.png

    Left: Stanford bunny without modifiers. Right: Stanford bunny with the "Fill holes" effect.
