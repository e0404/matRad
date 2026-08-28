matRad API Documentation
========================

matRad is separted into a top-level API that consists of a few functions to run the treatment planning workflow resting on the top-level of the the matRad folder / module. The core implementation is organized in several subdirectories / submodules that are called by the API functions and can be used for more fine-grained development.
Startup and configuration functions are located in the root of the repository.

.. note::

   The root-level scripts (``matRad.m``, ``matRadGUI.m``, ``matRad_rc.m`` etc.) are not auto-documented here.
   See the :ref:`run_script` guide for usage information.

Below, the top-level matRad functions are explained. For more specialized functions, refer to the documentation of the respective `Modules / Subfolders`_.

.. contents::
   :local:
   :depth: 1

Global Configuration
--------------------
.. _config:

matRad's global configuration class :class:`MatRad_Config` is used to set up the environment and configuration for the matRad application.
It is implemented as a Singleton pattern and thus consistent throughout a matRad session.

At the core, the class handles user folders, caches the environment (Matlab/Octave), provides matRad's version, and stores default parameters.
The class also provides logging functionality enabling control over output via log levels.

----

.. autoclass:: matRad.MatRad_Config
    :members:
    :undoc-members:
    :show-inheritance:
    :private-members:
    :noindex:
..

Top-level API functions
-----------------------

.. automodule:: matRad
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:
   :exclude-members: MatRad_Config

Modules / Subfolders
--------------------

.. toctree::
   :maxdepth: 5
   :glob:
   :includehidden:
   :reversed:

   *
   optimization/index

..
