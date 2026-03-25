.. |matRad_logo| image:: ../../matRad/gfx/matRad_logo.png
   :width: 80 px
   :alt: matRad
   :target: https://www.matRad.org

.. _quickstart:

Quick Start
===========

It's the first time you want to use matRad?

First, get a local copy of matRad by download or git cloning. Having done that, we recommend you navigate into the folder in Matlab and execute

.. code-block:: matlab

    matRad_rc

which will setup the path & configuration and tell you the current version.

Then there're three options for a pleasant start with matRad. Choose one or try out each of them.


.. rubric:: Option 1: Using the GUI
    :heading-level: 2


For an intuitive workflow with the graphical user interface, type

.. code-block:: matlab

    matRadGUI

in your command window. An empty GUI should be opened. Click the *Load.mat* data-Button in the Workflow-section to load a patient. Set the plan and optimization parameters, calculate the dose influence matrix and execute the fluence optimization in the GUI.

.. rubric:: Option 2: Using the main script
    :heading-level: 2

If you prefer scripting, open the default script *matRad.m* from the main matRad folder

.. code-block:: matlab

    edit matRad.m

Use it to learn something about the code structure and execute it section by section.

You can also run the full script for an example photon plan by just typing

.. code-block:: matlab

    matRad

in your command window.

.. rubric:: Option 3: Using the examples
    :heading-level: 2

The most time consuming but also most educational approach to matRad.

When in the main matRad folder, navigate to the folder *examples*. Open one of the examples given there. Execute it section by section. Move on to the next example afterwards.

**See also:**

.. toctree::
   :maxdepth: 1

   guiintro
   matradscript

..

