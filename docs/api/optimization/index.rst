############
optimization
############

.. toctree::
   :maxdepth: 4

   basefuncs
   optimizers
   projections

Optimization Problems
---------------------

An optimization problem combines the individual objectives into an objective function and organizes the constraint structure, handing a standard optimization problem to the :ref:`optimizer <optimizers>`.
Advanced planning problems are implemented as subclasses: direct aperture optimization for static apertures, and its VMAT counterpart for arcs.

.. autoclass:: matRad.optimization.matRad_OptimizationProblem
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:

.. autoclass:: matRad.optimization.matRad_OptimizationProblemDAO
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:

.. autoclass:: matRad.optimization.matRad_OptimizationProblemVMAT
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:

Local functions
---------------

.. automodule:: matRad.optimization
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:

..
