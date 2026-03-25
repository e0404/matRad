.. |matRad_logo| image:: ../../matRad/gfx/matRad_logo.png
   :width: 80 px
   :alt: matRad
   :target: https://www.matRad.org

.. _techdoc:

===============
Technical Guide
===============

|matRad_logo| features a very modular and sequential design which is reflected in the matRad script.
After importing your own data or loading one of the provided cases, you can start working with matRad dose calculation and optimization modules.
The four main parts of the matRad workflow are

.. image:: /images/matRad_steps.png

Information about the individual modules is given in the following sections:

.. toctree::
   :maxdepth: 2
   :hidden:

   config
   plan
   dosecalc
   visualization

Global configuration with :ref:`MatRad_Config <config>` 

:ref:`Set treatment plan parameters <plan>`

:ref:`Dose influence matrix calculation <dosecalc>`

:ref:`Fluence optimization <plan_opt>` (potentiall followed by sequencing)

:ref:`Visualization <visualization>`

Important variables and data structures
---------------------------------------

.. toctree::
   :maxdepth: 1
   :glob:

   ../datastructures/*
   ../datastructures/basedata/*

Additional information
----------------------

.. toctree::
   :maxdepth: 2
   :hidden:

   get
   requirements
   guioverview
   plan
   dosecalc
   optimization
   visualization
   coords
   dicomimport
   octave

.. toctree::
   :maxdepth: 2
   :hidden:

   ../algorithms/dicom
   ../algorithms/dosecalc
   ../algorithms/doseCalc/doseengines
   ../algorithms/optimization
   ../algorithms/sequencing



:ref:`How to run matRad with Octave <octave>`

:ref:`matRad coordinate system <coords>`

:ref:`The CORT dataset <cort>`

:ref:`DICOM import <dicomimport>`

:ref:`Minimum system requirements <requirements>`
