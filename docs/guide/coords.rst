.. _coords:

============================
The matRad Coordinate System
============================

.. _LPS:

LPS Coordinates
---------------

matRad uses the LPS (Left, Posterior, Superior) coordinate system.

.. image:: /images/CoordinateSystem/LPScoordinates.png
    :width: 400px

In this system, the x-axis points towards the left patient-side, the y-axis towards the posterior direction, and the z-axis towards the superior direction (see image below).

.. image:: /images/CoordinateSystem/LPS2.png
    :width: 400px

3D points, e.g., the isocenter in the ``pln`` struct or source and target points in the ``stf`` struct, directly follow the conventions of the LPS coordinate system: The first coordinate, e.g., ``pln.isoCenter(1)`` corresponds to the x coordinate (right-left direction). 
For the dose and CT cubes, which are stored as MATLAB 3D arrays, the second dimension corresponds to the x-coordinate (right-left direction), and the first dimension corresponds to the y-coordinate (anterior-posterior direction). 
This permutation is due to MATLAB's standard way of displaying two-dimensional matrices with the `image <http://de.mathworks.com/help/matlab/ref/image.html>`_ command, which displays the first array dimension along the vertical direction.

.. tip::

    In short this means that coordiantes (x/y/z) correspond to a (j,i,k) indexing.

    This is similar to the difference between MATLAB's `meshgrid <http://de.mathworks.com/help/matlab/ref/meshgrid.html>`_ and `ndgrid <http://de.mathworks.com/help/matlab/ref/ndgrid.html>`_ functions, where the first function uses the first dimension for the x-coordinate and the second dimension for the y-coordinate, while the second function uses the first dimension for the y-coordinate and the second dimension for the x-coordinate.

.. _voxelCoords:

World vs Cube System
---------------------

matRad differentiates between a world system and a cube system.

Coordinates in the :ref:`ct <ct>` struct as well as the plan isocenter, for example, are always given in world coordinates.

Cube coordinates are defined on an image cube (e.g. the :ref:`ct.cube <ct.cube>`). 
To enable fast access, they are defined such that the coordinates of the most right, anterior, inferior voxel center in the cube has the coordinates (``resolution.x`` / ``resolution.y`` / ``resolution.z``).
Consequently, the most right, anterior, inferior corner of the most right, anterior, inferior voxel is **not** located at (0 | 0 | 0) but at (resolution.x/2 | resolution.y/2 | resolution.z/2). 

.. danger::

    In versions prior to matRAd 3.1.0, the isocenter was given in the cube system. This was changed because during resampling of the dose grid, the isocenter often was different in the dose cubes coordinate system.

Conversion between world coordinates, cube coordinates and voxel indices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
matRad provides helper functions in the :mod:`geometry <matRad.geometry>` folder.

- :func:`matRad_world2cubeCoords` converts world coordinates to cube coordinates.
- :func:`matRad_cubeCoords2worldCoords` converts cube coordinates to world coordinates.
- :func:`matRad_world2cubeIndex` converts world coordinates to voxel indices.
- :func:`matRad_cubeIndex2worldCoords` converts voxel indices to world coordinates.

There is no function to convert cube coordinates to voxel indices, as this is a simple operation:

.. code-block:: matlab
    
    %Get the indices from coordinates
    coords = round(coords ./ [gridStruct.resolution.x gridStruct.resolution.y gridStruct.resolution.z]);

    %Do the permutation
    indices = coords(:,[2 1 3]);

.. _rotation:

Gantry and Couch Rotation
-------------------------

The rotation of the gantry (Φ) and the couch (θ) are defined as follows:

- **Gantry:** Clockwise rotation around the z-axis
- **Couch:** Counter-clockwise rotation around the y-axis

The rotation matrix can be obtained with :func:`matRad_getRotationMatrix`.

Simple gantry rotation:

.. image:: /images/CoordinateSystem/RotatingGantry.gif
    :width: 400px

Simple couch rotation:

.. image:: /images/CoordinateSystem/RotatingCouch.gif
    :width: 400px

Simultaneous couch and gantry rotation:

.. image:: /images/CoordinateSystem/RotatingCouch+Gantry.gif
    :width: 400px
