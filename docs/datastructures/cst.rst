.. _cst:

==================
The cst cell Array
==================

The constraints of all defined volumes of interest (VOIs) are stored inside the ``cst`` cell array. It is structured as follows:

.. _cst-cell:

Screenshot of the cst-cell:
.. image:: /images/cstCellScreenshot.png
    :alt: Screenshot of the cst cell

.. list-table:: Structure of the cst cell array
    :header-rows: 1

    * - Column
      - Content
      - Description
    * - **1**
      - :ref:`VOI index <VolInd>`
      - Number to identify the VOI
    * - **2**
      - :ref:`VOI name <VolName>`
      - String describing the VOI
    * - **3**
      - :ref:`VOI type <VolType>`
      - Specification whether the VOI is an organ at risk (OAR), a target volume or should be ignored
    * - **4**
      - :ref:`Voxel indices <VoxInd>`
      - Vectors containing the indices of all voxels of the CT that are covered by the VOI. Stored as a cell array of vectors (for enabling handling of multiple scenarios)
    * - **5**
      - :ref:`Tissue parameters <TissParam>`
      - Structure containing information about the tissue of the VOI and its overlap priority
    * - **6**
      - :ref:`Dose objectives <DoseParam>`
      - Cell array containing information about the functions used to calculate the objective & constraint function value
    * - **7**
      - Precomputed Contours
      - After GUI startup, this column contains precomputed contour data for display

.. _VolInd:

VOI index
---------

All defined VOIs are enumerated starting with 0.

.. _VolName:

VOI name
--------

The VOI name is a string containing an organ name or a short description of the volume (e.g. ``BODY``, ``Liver``, ``GTV``, ...).

.. _VolType:

VOI type
--------

The VOI type specifies how the volume is considered during treatment planning:

.. list-table:: VOI types and their handling during treatment planning
    :header-rows: 1

    * - VOI type
      - Handling during treatment planning
    * - **TARGET**
      - The VOI will be covered with spot positions (protons / carbon ions) and bixels (photons) as defined in the :ref:`stf struct <stf>`. During the fluence optimization, it will be considered according to the defined :ref:`dose objectives <DoseParam>`.
    * - **OAR**
      - The VOI will not be covered with spot positions or bixels. During the fluence optimization, it will be considered according to the defined :ref:`dose objectives <DoseParam>`.
    * - **IGNORED**
      - The VOI will not be considered during the treatment planning.

.. _VoxInd:

Voxel indices
-------------

The indices of all voxels (of the :ref:`CT-cube <ct>`) that are covered by the VOI are stored in a vector within a cell array. I.e. we store the segmentation for the VOI as a binary mask, the polygon contour data is not part of matRad's standard data sets.
As the same voxel can be covered by more than one VOI, an overlap priority (see :ref:`tissue parameters <TissParam>`) is defined to handle potential discrepancies when calculating the objective function value and generating the :ref:`stf struct <stf>`.

.. _TissParam:

Tissue parameters
-----------------

.. image:: /images/cstCellTissueParametersScreenshot.png
    :alt: Screenshot of tissue parameters

Data can also be stored as in the :ref:`old format (see below) <DoseParamOld>`.

New constraints or objectives can be implemented by adding a respective class definition to the :mod:`matRad.optimization.+DoseConstraints` or :mod:`matRad.optimization.+DoseObjectives` folder.

.. _DoseParam:

Dose Objectives & Constraints since v2.10.0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
matRad supports inverse planning based on the minimization of a weighted sum of objectives subject to non-linear yet differentiable hard constraints. The following kind of individual objectives are currently supported:

.. include:: ../includes/objtable.rst

Constraints are somewhat built around similar goals as obejctives:

.. include:: ../includes/constrtable.rst

When generating an objective / constraint from code, we suggest to wrap the instantiation of the objective/constraint in a ``struct()`` call, as shown in the first phantom example:

.. literalinclude:: ../../examples/matRad_example1_phantom.m
    :caption: examples/matRad_example1_phantom.m
    :lines: 47-49
    :lineno-match:
    :language: matlab


This will ensure that, when saving to a mat-file, we don't save the class object, which improves compatibility.

.. _DoseParamOld:

Before Version 2.10.0
~~~~~~~~~~~~~~~~~~~~~
In the earlier version, matRad stored the objectives and constraints defined for inverse planning as an array of structs. :func:`matRad_convertOldCstToNewCstObjectives` can be used to convert the old definition to the new format.

.. _defaultValues:

Default *cst*-values
--------------------

The patient data contained within matRad (*ALDERSON, BOXPHANTOM, HEAD_AND_NECK, LIVER, PROSTATE and TG119*) have default values defined within the *cst*-cell. 
 
These values are chosen to produce a reasonable treatment plan, when using coplanar and equidistant photon beams. They can be used as a reference point for more sophisticated treatment plans.