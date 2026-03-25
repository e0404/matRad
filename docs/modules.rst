matRad API Documentation
========================

matRad is separted into a top-level API that consists of a few functions to run the treatment planning workflow resting on the top-level of the the matRad folder / module. The core implementation is organized in several subdirectories / submodules that are called by the API functions and can be used for more fine-grained development.
Startup and configuration functions are located in the root of the repository.

.. note::

   The root-level scripts (``matRad.m``, ``matRadGUI.m``, ``matRad_rc.m`` etc.) are not auto-documented here.
   See the :ref:`run_script` guide for usage information.

Below, the top-level matRad functions are explained. For more specialized functions, refer to the documentation of the respective `Modules / Subfolders`_.

Top-level API functions
-----------------------

.. automodule:: matRad
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:

.. toctree::
   :maxdepth: 3
   api/config
   
Modules / Subfolders
--------------------

.. toctree::
   :maxdepth: 3
   :glob:
   :includehidden:
   
   api/**

..
