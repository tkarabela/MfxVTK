Sample Points (Volume)
**********************

Generate random points inside given mesh.

.. tip::
    In Blender, parent an object to the resulting point cloud and set
    "Instancing" to "Verts" in the point cloud Object properties.

    This way you can create similar effects to a particle system with
    *Emit from volume*. Points from the MfxVTK effect have a different look
    and they are static -- you don't have to worry about generating particles
    or baking them.

:Input: polygonal mesh (should be closed)
:Output: point cloud
:VTK classes: ``vtkImplicitPolyDataDistance``

Options
#######

.. figure:: /_static/properties-sample-points-volume.png
    :align: right
    :width: 300px

Number of points
    How many points should be generated.

    .. note::
       If your mesh is very thin or tiny relative to its bounding box,
       the effect may fail to sample all points. In this case, increasing
       desired number of points may help (it will not give up as early).

Distribute uniformly
   If this option is turned on, the points will be distributed pretty evenly across
   the volume (quasi-random, low discrepancy sampling).

   If this option is off, the points will be sampled more chaotically
   across the volume (pseudo-random uniform distribution over bounding box).

Auto simplify
   Allow the effect to simplify (decimate) the mesh for sampling, dramatically
   improving performance for large meshes. In rare cases, this may lead to points
   being sampled outside the original mesh. If this happens, you can try turning
   this off.

Example
#######

.. figure:: /_static/example-sample-points-volume.png

    Left: The Stanford bunny.
    Middle: "Sample points (volume)" effect with "Distribute uniformly" on.
    Right: "Sample points (volume)" effect with "Distribute uniformly" off.
