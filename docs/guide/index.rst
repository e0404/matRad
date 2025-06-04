.. |matRad_logo| image:: ../../matRad/gfx/matRad_logo.png
   :width: 80 px
   :alt: matRad
   :target: https://www.matRad.org

.. _techdoc:

=======================
Technical Documentation
=======================

|matRad_logo| features a very modular and sequential design which is reflected in the matRad script. 
After importing your own data or loading one of the provided cases, you can start working with matRad dose calculation and optimization modules. 
The four main parts of the matRad workflow are

.. image:: /images/matRad_steps.png

Information about the individual modules is given in the following subsections:

:ref:`Set treatment plan parameters <plan>`

:ref:`Dose influence matrix calculation <dosecalc>`

:ref:`Fluence optimization <plan_opt>` (potentiall followed by sequencing)

:ref:`Visualization <visualization>`

How to cite matRad
------------------

:ref:`matRad publications <cite>`

.. _httpsrawgitcomwikie0404matradimagesmatrad_blanksvg--height--25pxs-most-important-matlab-variables:

.. _matrad_variables:

|matRad_logo|'s most important Matlab variables
-----------------------------------------------

:ref:`pln-struct <pln>` Treatment plan information

:ref:`ct-struct <ct>` CT-data

:ref:`cst-cell array <cst>` Structure sets, inverse planning objectives, and other meta-information

:ref:`stf-struct <stf>` Steering information

:ref:`dij-struct <dij>` Dose influence data

:ref:`result-struct <result>` Resulting dose distribution, RBE cube etc..

.. toctree::
   :maxdepth: 2
   :hidden:
   :glob:

   datastructures/*

|matRad_logo|'s top-level functions
-----------------------------------

:ref:`MatRad_Config <config>` Global configuration class

.. toctree::
   :maxdepth: 2
   :hidden:
   :glob:

   api/*

Additional information
----------------------

:ref:`How to run matRad with Octave <octave>`

:ref:`matRad coordinate system <coords>`

:ref:`The CORT dataset <cort>`

:ref:`DICOM import <dicomimport>`

:ref:`Minimum system requirements <requirements>`

