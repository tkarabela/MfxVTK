Poke
****

This effect performs softbody-like deformation, displacing input mesh to avoid
collision with collider mesh. Since it has no internal physics or time dependence,
it can be easily posed/animated. It works with one or two input meshes.

.. tip::
   You can use the color attribute to specify which parts of the input and/or
   collider mesh should be taken into account. When a color attribute is present,
   only **white** (RGB = 255, 255, 255) **areas of the mesh** and only **black**
   (RGB = 0, 0, 0) **areas of the collider** participate in the collision.

.. note::
   If you **don't specify collider mesh**, the color attribute is **mandatory** (so that
   the effect knows which parts are supposed to be the collider ie. "rigid", and which
   parts should be deformed ie. "soft".)

   To prevent artifacts, it's important to keep "rigid" and "soft" areas from touching
   along the surface; there should not be any polygons with both "rigid" and "soft" vertices. Leave a
   gray area between the white and black parts.

:Input: 1 or 2 polygonal meshes (with optional color attribute)
:Output: polygonal mesh
:VTK classes: ``vtkStaticCellLinks``, ``vtkModifiedBSPTree``, ``vtkPolyDataNormals``
:Preserves topology: Yes
:Multithreaded: Yes (only building topology)

Options
#######

.. figure:: /_static/properties-poke.png
    :align: right
    :width: 300px


Maximum distance
    Determines how far away from collider can a collision be detected.
    Default value 0.0 means auto (bounding box diagonal). Small values
    may cause the collider to sink into mesh and not detect the collision.
    Large values increase computation time.

Falloff radius
    Controls area around deformation which should be smoothed. Large values
    create larger area of depression around points of contact.

Falloff exponent
    Controls area around deformation which should be smoothed. Small values
    create sharper edge at points of contact, while larger values make the
    transition more smooth.

Collision smoothing ratio
    Controls smoothing inside depressed area. Larger values make the displaced
    surface more smooth, but may reintroduce collision as the mesh is stretching
    back.

Collider normal factor
    Blends between displacement along mesh and collider normals.

Number of iterations
    Controls Laplacian smoothing strength; more iterations make smoother surface,
    may reintroduce collisions as the mesh is stretching back, and increase
    computation time.

Offset
    Make the collider larger for collision detection (ie. like if it was magnetic,
    repulsing the mesh at a distance). This can be used to counteract smoothing
    when it reintroduces collisions.

Debug
    (For development) Save mesh and collider in VTK format for debugging.

.. figure:: /_static/poke-falloff-exponent.png

    Effect of "Falloff exponent"

.. figure:: /_static/poke-falloff-radius.png

    Effect of "Falloff radius"

.. figure:: /_static/poke-iterations.png

    Effect of "Number of iterations"

.. figure:: /_static/poke-offset.png

    Effect of "Offset"

.. figure:: /_static/poke-collision-smoothing.png

    Effect of "Collision smoothing ratio"

Example
#######

.. figure:: /_static/example-poke.png

    Suzanne laying on a Plane with applied "Poke" filter.
