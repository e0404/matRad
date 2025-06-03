.. _config:

The global configuration class
==============================

matRad's global configuration class is used to set up the environment and configuration for the matRad application. 
It is implemented as a Singleton pattern and thus consistent throughout a matRad session.

At the core, the class handles user folders, caches the environment (Matlab/Octave), provides matRad's version, and stores default parameters. 
The class also provides logging functionality enabling control over output via log levels.

----

.. autoclass:: matRad.MatRad_Config
    :members:
    :undoc-members:
    :show-inheritance:
    :private-members:
