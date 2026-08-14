.. _obj_constr:

======================
Optimization Functions
======================

matRad implements objectives and constraints for the optimization of dose distributions, as well as objectives acting directly on the fluence.
All of them derive from :class:`matRad_OptimizationFunction`, which provides what is independent of the optimized quantity: the parameter description used by the GUI, struct serialization and instantiation from a struct.

From it derive the quantity specific bases: :class:`matRad_DoseOptimizationFunction` for functions of the dose, differentiating between dose :ref:`objectives <objectives>` and :ref:`constraints <constraints>`, and :class:`matRad_FluenceOptimizationFunction` for :ref:`fluence objectives <fluence_objectives>`.

.. _objectives:

Dose Objectives
---------------

.. automodule:: matRad.optimization.+DoseObjectives
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:

.. _constraints:

Dose Constraints
----------------

.. automodule:: matRad.optimization.+DoseConstraints
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:

.. _fluence_objectives:

Fluence Objectives
------------------

Fluence objectives act on the bixel weight vector directly instead of on a dose related quantity, so that they take no part in the backprojection.
They are given either in ``pln.propOpt.fluenceObjectives`` or alongside the dose objectives in the :ref:`cst cell <cst>`, in which case the structure they are attached to is ignored - the objective still acts on the whole fluence.
Since they only need the bixel to beam mapping from the dose influence matrix, they work for photons and particles alike.

.. automodule:: matRad.optimization.+FluenceObjectives
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:
