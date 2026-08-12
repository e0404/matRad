.. _dosecalc:

########
doseCalc
########

.. contents::
   :local:

Dose calculation engines and supporting utilities for photon, proton, and ion therapy.

.. automodule:: matRad.doseCalc
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:

Dose Engines
------------

.. automodule:: matRad.doseCalc.+DoseEngines
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:

The base class of the engine hierarchy and the FRED engine live in class folders, which the directive above does not descend into:

.. autoclass:: matRad.doseCalc.+DoseEngines.matRad_DoseEngineBase
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:

.. autoclass:: matRad.doseCalc.+DoseEngines.matRad_ParticleFREDEngine
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:

FRED
----

Interface components for the FRED Monte Carlo dose calculation engine.

.. automodule:: matRad.doseCalc.FRED
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:

MCsquare
--------

Interface components for the MCsquare fast Monte Carlo proton dose calculation engine.

.. automodule:: matRad.doseCalc.MCsquare
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:

TOPAS
-----

Interface componentsfor the TOPAS Monte Carlo dose calculation engine.

.. automodule:: matRad.doseCalc.topas
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:

..
