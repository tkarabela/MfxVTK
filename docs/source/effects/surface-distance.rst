Distance Along Surface
**********************

This effect computes shortest distance along surface to given area(s) on the mesh.
You can imagine this as heating some part of the mesh to glowing red temperature,
this effect computes how the temperature gradually changes along the surface into
the cooler areas.

To use this effect, create a vertex color on your mesh and set all vertices to black.
Then paint some vertices white -- these will be the starting ("glowing hot") points.
It's okay to use multiple points or even multiple "islands" on different parts of the mesh.

Output of the effect is a UV map, where the U component contains the distance.

.. tip::
    You can then use the UV map in your shader with a color ramp to create
    the desired effect (a wave running along the surface, mesh gradually
    appearing, etc.).

.. note::
    This is an approximation (upper bound) of the true distance -- for best results,
    use quads and even topology. Polygons with a high aspect ratio or
    sudden changes in point density may cause artifacts (the distance will appear longer than it really is).

:Input: polygonal mesh (with color attribute)
:Output: polygonal mesh (with UV map)
:VTK classes: ``vtkStaticCellLinks``
:Preserves topology: Yes
:Multithreaded: Yes (only building topology)

Options
#######

.. figure:: /_static/properties-surface-distance.png
    :align: right
    :width: 300px


Normalize distance
    If this option turned on, computed distance in the U component will be
    scaled to [0;1] range. This is a safe default, suitable for cases when
    you don't need to match distances across multiple objects.

    If this option is turned off, the U component will be set to the actual
    distance along surface (as computed from local model coordinates). This makes
    distances comparable across objects. Note that this will generally cause
    the U component to have values greater than 1.0,
    which may cause problems in your 3D application.


Example
#######

.. figure:: /_static/example-surface-distance.png

    Left: Stanford bunny with input color attribute (source vertices are white,
    rest of the mesh is black).
    Right: Stanford bunny shaded with U component from *Distance Along Surface* UV map
    (red is 0.0, blue is 1.0).
