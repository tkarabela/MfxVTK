List of Effects
---------------

On this page, you can find reference to all available effects in MfxVTK.
They are roughly divided into three categories:

- effects related to **point clouds** (meshes with no faces or edges),
- effects related to **edge wireframes** (meshes with no faces) and
- effects related to **polygonal meshes** (the usual kind of mesh).

Depending on the effect, different kinds of input are supported
(eg. :doc:`effects/make-tubes` can work with wireframe *or* polygonal mesh, while
:doc:`effects/fill-holes` *only* works with polygonal mesh). If a given effect
doesn't support the input you give it, it will either have no effect or
you may get empty output. Effects should generally "just work" for all combinations
that make sense.

.. note::
   It's totally OK and even encouraged to **stack the effects** on top of one another.
   You may start with a polygonal mesh, create a point cloud out of it, then create
   a wireframe, then a polygonal mesh out of the wireframe.

   The possibilities are endless! (Host limitations notwithstanding; Blender lends
   itself to work with polygonal meshes more so than with point clouds or wireframes,
   though all are possible.)

.. toctree::
   :maxdepth: 1
   :caption: Point Cloud Effects

   effects/sample-points-surface
   effects/sample-points-volume

.. toctree::
   :maxdepth: 1
   :caption: Edge Wireframe Effects

   effects/extract-edges
   effects/make-tubes
   effects/tetrahedral-wireframe

.. toctree::
   :maxdepth: 1
   :caption: Polygonal Mesh Effects

   effects/smooth
   effects/decimate
   effects/fill-holes

.. Connect Tetrahedra
    ******************
    vtkDelaunay3D

.. Decimate Polyline
    *****************
    decimate polyline

.. Build Tree
    **********
    spanning tree from points/wireframe


.. Vertex Group to UV
    ******************
    Blender-specific effect, put vertex group to UV.


.. Reconstruct Surface
    *******************
