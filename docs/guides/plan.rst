.. _plan:

===========================
How to configure your plan?
===========================

Before you can start with dose calculation and inverse planning in matRad you need to set a couple of general treatment plan parameters.

Within the main matRad script, this corresponds to the first cell. Within the GUI, these settings can be adjusted interactively.

.. literalinclude:: ../../matRad.m
    :lines: 18-53
    :lineno-match:
    :language: matlab

After import of patient data in matRad's native format, the desired settings are specified in the :ref:`pln struct <pln>`. 

Note that besides some basic parameters (like number of fractions or the machine), there are workflow-specific property structures in ``pln.prop*``. These structures contain the parameters for the different workflows, such as the dose influence matrix calculation, optimization, and sequencing. matRad will instatiate the appropriate algorithm based on these configurations and try to override their properties.

This script uses the top-level API exposed when importing pyRadPlan. The top-level functions are designed to take the main data structures as input and configure the corresponding workflow step via the `Plan` using the attribute dictionaries `pln.prop*`:

.. include:: ../includes/planapi.rst

Please see the corresponding page about the :ref:`pln struct <pln>` for further information.

Note that matRad comes with sample patient data from the `CORT dataset: common optimization for radiation therapy <https://academic.oup.com/gigascience/article/3/1/2047-217X-3-37/2682969>`_ so you can directly play around with all functionalities. More information can be found on the  :ref:`dedicated CORT dataset section <cort>`.