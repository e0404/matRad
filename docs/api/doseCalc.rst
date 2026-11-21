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

MCNP
----

Interface components and data for the MCNP Monte Carlo neutron dose calculation engine
(see :class:`DoseEngines.matRad_NeutronMCNPEngine`). The helper functions in this module and
its ``tools`` submodule are private to the engine and therefore not part of the API reference.

.. automodule:: matRad.doseCalc.MCNP

.. autoclass:: matRad.doseCalc.MCNP.matRad_MCNPConfig
   :members:
   :undoc-members:
   :show-inheritance:

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
