Example Cookbook
================

This page shows examples of using MfxVTK in Blender.

Tetrahedral Bunny
-----------------

.. figure:: /_static/bunny-tetrahedral-mesh.png
    :align: center

- *Sample points (volume)*
   - 3000 points
   - *Distribute uniformly* turned on
- *Tetrahedral wireframe*
   - Maximum edge length = 0.4
- *Make tubes*
   - Radius = 0.02
   - Number of sides = 6

Bunny from Spheres
------------------

.. figure:: /_static/bunny-spheres.png
    :align: center

- *Sample points (volume)*
   - 6000 points
   - *Distribute uniformly* turned off
- Create a UV sphere, parent it to Bunny, set *Instancing* to *Verts* in Bunny's Object properties
