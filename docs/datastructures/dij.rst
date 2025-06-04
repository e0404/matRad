.. _dij:

======================
The dij Data Structure
======================

The dij struct holds pre-calculated dose influence data for inverse planning. The individual fields contain the following information.

Starting from version 3.1.0, the dose grid has been separated from the CT grid, meaning that each can have a separate resolution depending on the application. Therefore, the two substructures ``doseGrid`` and ``ctGrid`` hold the following information for both the CT and the dose grid.

**dij.doseGrid.resolution - dij.ctGrid.resolution**

    The resolution of an individual voxel in the dose/ct cube in [mm] in x-, y-, and z-direction (Note that the default value for the doseGrid resolution is 2.5 mm, 2.5 mm and 3 mm in x-, y- and z- directions, and is not assigned equal to the ct resolution).

**dij.doseGrid.dimensions - dij.ctGrid.dimensions**

    The resulting dose/ct cube dimension in voxels.

**dij.doseGrid.numOfVoxels - dij.ctGrid.numOfVoxels**

    The total number of voxels in the entire dose/ct cube.

**dij.numOfBeams**

    Specifies the number of beams used for dose calculation and consequently inverse planning.

**dij.numOfScenarios**

    Number of scenarios considered during dose calculation. Usually, only one scenario is used. In a special research mode, however, it is possible to compute dose on multiple 4D CT phases, considering isocenter shifts or range uncertainties.

**dij.numOfRaysPerBeam**

    Specifies the number of rays per beam. For photons, this number also corresponds to the number of bixels per beam. For particles however, it is also possible to have multiple spot positions with different energies at the same lateral spot position.

**dij.totalNumOfRays**

    Specifies the total number of all rays, i.e. ``dij.totalNumOfRays = sum(dij.numOfRaysPerBeam)``.

**dij.totalNumOfBixels**

    Specifies the total number of bixels over all beams.

**dij.bixelNum**

    Lists the bixel number in an individual beam for all columns in the precomputed influence data. Together with ``dij.rayNum`` and ``dij.beamNum`` this information facilitates an easy assignment of columns of the influence data to the :ref:`stf struct`.

**dij.rayNum**

    Lists the ray number in an individual beam for all columns in the precomputed influence data. Together with ``dij.bixelNum`` and ``dij.beamNum`` this information facilitates an easy assignment of columns of the influence data to the :ref:`stf struct`.

**dij.beamNum**

    Lists the beam number for all columns in the precomputed influence data. Together with ``dij.bixelNum`` and ``dij.rayNum`` this information facilitates an easy assignment of columns of the influence data to the :ref:`stf struct`.

**dij.physicalDose**

    Pre-computed dose influence matrix with ``dij.doseGrid.numOfVoxels`` rows and ``dij.totalNumOfBixels`` columns stored within a cell array (Note that starting from Version 3.0, this matrix is built by default based on ``numOfVoxels`` on doseGrid, not the ctGrid). This matrix specifies the dose contribution from every bixel to every voxel, stored with MATLAB's built-in double precision sparse matrix format.

**dij.mAlphaDose**

    Pre-computed alpha*dose matrix with ``dij.numOfVoxels`` rows and ``dij.totalNumOfBixels`` columns, stored with MATLAB's built-in double precision sparse matrix format within a cell array. This matrix is only computed for biological optimization, where this information is required to compute dose-averaged alpha cubes that are in turn required for three-dimensional RBE modeling.

**dij.mSqrtBetaDose**

    Pre-computed sqrt(beta)*dose matrix with ``dij.numOfVoxels`` rows and ``dij.totalNumOfBixels`` columns, stored with MATLAB's built-in double precision sparse matrix format within a cell array. This matrix is only computed for biological optimization, where this information is required to compute dose-averaged alpha cubes that are in turn required for three-dimensional RBE modeling.

Screenshot of the dij-struct:
    
.. image:: /images/dij-struct.png
   :alt: Screenshot of the dij struct
