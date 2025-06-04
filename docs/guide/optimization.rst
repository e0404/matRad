.. _plan_opt:

====================
Fluence Optimization
====================

The goal of the fluence optimization is to find a set of bixel/spot weights that yield the best possible dose distribution according to the clinical objectives and constraints underlying the radiation treatment.

For mathematical optimization, these clinical objectives and constraints have to be translated into mathematical objectives and constraints. matRad supports the mathematical optimization of a weighted sum of objectives to help finding an optimal trade-off between adequate target coverage and normal tissue sparing for an individual patient as well as the formulation of constraints. The individual objectives and constraints are defined per structure, can be chosen by the user.

The overall fluence optimization process is coordinated by the top-level function :func:`matRad_fluenceOptimization`.

Since Version 2.10.0
--------------------

*Since version 2.10.0 objectives and constraints are implemented with an object oriented approach. For the old format, see further down below.*

At the moment, matRad allows for objectives and constraints based on dose. Each objective and constraint is defined in a class derived from :class:`matRad_DoseOptimizationFunction`. Objectives and Constraints are distinguished by the abstract subclasses :class:`DoseObjectives.matRad_DoseObjective` & :class:`DoseConstraints.matRad_DoseConstraint` within the respective package folders :mod:`. This enables the easy implementation of new (dose-based) constraints & objectives by writing your own class, inheriting from those and implementing the therein declared interface, i.e. the parameter definition, and respective objective/constraint function and its gradient/jacobian.
New objectives/constraints are then automatically recognized also in the GUI. 

Currently, the available objectives are found in the :mod:`matRad.optimization.+DoseObjectives` package and the available constraints in the :mod:`matRad.optimization.+DoseConstraints` package. The available objectives and constraints are listed in the following tables:

.. include:: ../includes/objtable.rst

.. include:: ../includes/constrtable.rst

This is extended to the implementation of optimization problems and :mod:`optimizers <matRad.optimization.optimizer>`. An optimization problem :class:`matRad_OptimizationProblem` combines the single objectives into an objective function and organizes the constraint structure to give a standard optimization problem to the optimizer. matRad implements optimizers as derived classes of :class:`matRad_Optimizer` and defaults to the `IPOPT <https://www.coin-or.org/Ipopt/>`_ package for large scale non-linear optimization which is included via a MEX file and interfaced in :class:`matRad_OptimizerIPOPT`. MATLAB's proprietary fmincon support is added in :class:`matRad_OptimizerFmincon`, other optimizers may be added as well. So far, only non-linear constrained optimization is supported by :class:`matRad_OptimizationProblem` and for optimizers. 
Optimizers can be changed by setting ``pln.propOpt.optimizer``.
The :class:`matRad_OptimizationProblem` class also enables to implement advanced planning problems as subclasses, like direct aperture optimization as implemented in :class:`matRad_OptimizationProblemDAO`.

Before Version 2.10.0
---------------------

The objectives and constraints are stored as a :ref:`dose objective struct <DoseParam>` within the :ref:`cst cell <cst>`. The objectives and constraints can be set including all necessary parameters via the :ref:`matRad GUI <run_gui>`. All functions involved in the optimization process are located in a subfolder called "optimization" within the matRad root folder. matRad relies on the `IPOPT <https://www.coin-or.org/Ipopt/>`_ package for large scale non-linear optimization which is included via a MEX file. IPOPT requires call back functions for objective function, gradient, constraint, and Jacobian evaluation. We use the wrapper functions ``matRad_objFuncWrapper.m``, ``matRad_gradFuncWrapper.m``, ``matRad_constFuncWrapper.m``, and ``matRad_jacobFuncWrapper.m`` to coordinate the evaluation of all defined objectives and constraints.

Optimization based on dose, effect, and RBE weighted dose
---------------------------------------------------------

All optimization functionalities work equally for optimization processes based on physical dose as well as biological effect `Wilkens & Oelfke (2006) <http://iopscience.iop.org/0031-9155/51/12/009>`_ and RBE-weighted dose according to `Krämer & Scholz (2006) <http://iopscience.iop.org/0031-9155/51/8/001>`_. The biological effect and the RBE-weighted dose are calculated with α and β base data that has been calculated according to the local effect model IV. α and β tables are available as part of the base data set `carbon_Generic <https://github.com/e0404/matRad/blob/master/matRad/basedata/carbon_Generic.mat>`_ which is provided with the matRad release.

Direct aperture optimization
----------------------------

For photons, matRad also features an experimental direct aperture optimization that largely follows the implementation described in `Wild et al. (2015) <https://doi.org/10.1118/1.4914863>`_ which is based on `Bzdusek et al. (2009) <http://www.ncbi.nlm.nih.gov/pubmed/19610322>`_ and (with some modification) `Unkelbach & Cassioli (2012) <http://iopscience.iop.org/article/10.1088/0031-9155/58/2/301/meta;jsessionid=B918D167B0BD2950E5A559F03F6CC517.c2.iopscience.cld.iop.org>`_.


