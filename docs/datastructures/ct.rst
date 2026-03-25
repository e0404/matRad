.. _ct:

=====================
The ct Data Structure
=====================

The ct structure contains, among others, the 3D-CT-cube (see `ct.cube`_), obtained from the planning-CT, and the voxel resolution (see `ct.resolution`_).

Screenshot of the ct-struct:
  
.. image:: /images/ctDataScreenshot.png

.. _ct-cube:

ct.cube
-------

The cube is a N\ :sub:`x` × N\ :sub:`y` × N\ :sub:`z` matrix (N\ :sub:`x,y,z` = number of voxels in x-, y- and z-direction) containing the water equivalent thickness of each voxel. We already translate HU to water equivalent electron density according to a look up table upon patient data import. The cube(s) is/are stored within a cell array to support multiple CT phases for 4D data.

.. _resolution:

ct.resolution
-------------

The resolution specifies the size of each voxel in x-, y-, and z-direction in [mm].

.. _cubeDim:

ct.cubeDim
----------

Number of voxels in x-, y- and z-direction (N\ :sub:`x`, N\ :sub:`y` and N\ :sub:`z`).

.. _numOfCtScen:

ct.numOfCtScen
--------------

The number of considered CT scenarios. Usually, this corresponds to one but in a special research mode, it is also possible to handle multiple CTs.