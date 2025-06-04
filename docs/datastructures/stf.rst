.. _stf:

======================
The stf Data Structure
======================

The stf struct holds all geometric information about the irradiation. The individual fields contain the following information.

Screenshot of the stf structure:
  
.. image:: /images/stfStructScreenshot.png
   :alt: stf structure screenshot

.. list-table::
    :header-rows: 1

    * - **Field**
      - **Description**
    * - **gantry angle**
      - Specification of the gantry angle for each beam in °. Range: 0° - 359°
    * - **couch angle**
      - Specification of the couch angle for each beam in °. Range: 0°- 359°
    * - **bixel width**
      - *Photons*: lateral width of each bixel. *Particles*: lateral spot spacing.
    * - **radiation mode**
      - Specification of the radiation mode. Options are: 'photons', 'protons' or 'carbon'.
    * - **number of Rays**
      - Number of rays for each beam
    * - **ray**
      - The :ref:`ray substructure <ray>` contains the information about the position and orientation of each ray / bixel for a single beam.
    * - **source point**
      - Position of the virtual radiation source in :ref:`LPS coordinates <matRad-coordinate-system>`.
    * - **number of bixels per ray**
      - *Photons*: each ray is a single bixel. *Particles*: :ref:`number of bixels / dose spots along each ray <numOfBixels>`.
    * - **total number of bixels**
      - Total number of bixels for each beam.

.. _ray:

stf.ray substructure
====================

The *stf.ray* substructure contains the information about the position and orientation of each ray/bixel for a single beam.

Screenshot of the stf.ray substructure:
  
.. image:: /images/stfStructRayScreenshot.png
   :alt: stf.ray substructure screenshot

.. list-table::
    :header-rows: 1

    * - **Field**
      - **Code**
      - **Description**
    * - **ray position** *(beam's eye view)*
      - ``rayPos_bev``
      - Point where the ray crosses the isocenter plane in the beam's eye view.
    * - **target point** *(beam's eye view)*
      - ``targetPoint_bev``
      - Target point of the ray from the beam's eye view. The target point extends the ray from source to rayPos_bev to behind the patient for use in rayTracing.
    * - **ray position** *(LPS-coordinates)*
      - ``rayPos``
      - rayPos_bev in the :ref:`LPS coordinate system <matRad-coordinate-system>`.
    * - **target point** *(LPS-coordinates)*
      - ``targetPoint``
      - targetPoint_bev in the :ref:`LPS coordinate system <matRad-coordinate-system>`.
    * - :ref:`energy <energy>`
      - ``energy``
      - *Photons*: Single value (max. LINAC energy). *Particles*: Energy values of the dose spots along the ray. The number of energies corresponds to the :ref:`number of bixels per ray <numOfBixels>`.

.. _energy:

stf.ray.energy field
====================

The *stf.ray.energy* field contains all energy values for the spots along the specified ray. If the radiation mode *'Photons'* was selected, it will be set to *NaN*.

Screenshot of the stf.ray.energy field:
  
.. image:: /images/stfStructRayEnergyScreenshot.png
   :alt: stf.ray.energy field screenshot

.. _numOfBixels:

stf.numOfBixelsPerRay field
===========================

The *stf.numOfBixelsPerRay* field contains the number of bixels / dose spots for every ray of the specified beam. This corresponds to the number of energy values per ray in the :ref:`stf.ray <ray>` substructure.

Screenshot of the stf.numOfBixelsPerRay field:
  
.. image:: /images/stfStructNumOfBixelsScreenshot.png
   :alt: stf.numOfBixelsPerRay field screenshot