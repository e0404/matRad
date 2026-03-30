.. _gui:

###
gui
###

.. contents::
   :local:


Main GUI
--------

The main graphical user interface of matRad is implemented in the class :class:`matRad_MainGUI`.
It provides access to all functionalities of matRad and is the central hub for visualization and interaction with the treatment planning workflow.

It is composed of several widgets that are implemented in the :mod:`matRad.gui.widgets` submodule.
The GUI can also be started with :func:`matRadGUI` from the root folder.

.. autoclass:: matRad.gui.matRad_MainGUI
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:

Basic GUI functions
-------------------

.. automodule:: matRad.gui
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:
   :exclude-members: matRad_MainGUI

Widgets
-------

.. automodule:: matRad.gui.widgets
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:
