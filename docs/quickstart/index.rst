.. |matRad_logo| image:: ../../matRad/gfx/matRad_logo.png
   :width: 80 px
   :alt: matRad
   :target: https://www.matRad.org

.. _quickstart:

Quick Start
===========

Is it the first time you want to use matRad?

First, get a local copy of matRad by downloading or git cloning. Having done that, we recommend you navigate into the folder in Matlab and execute

.. code-block:: matlab

    matRad_rc

which will setup the path & configuration and tell you the current version.

A key feature of matRad is that it can be used either through a graphical user interface (GUI) or as a MATLAB script.

Accordingly, there are different ways to familiarize yourself with the software and get started with matRad. Below, we present three different beginner-friendly approaches. Choose one or explore and try out each of them.

.. rubric:: Option 1: Using the GUI
    :heading-level: 2

For the most intuitive workflow, we recommend using the graphical user interface (GUI). In order to start it, type

.. code-block:: matlab

    matRadGUI

in your command window. An empty GUI should be opened. Click the *Load.mat* data-Button in the Workflow-section to load a patient. Set the plan and optimization parameters, calculate the dose influence matrix and execute the fluence optimization in the GUI.

.. rubric:: Option 2: Using the main script
    :heading-level: 2

If you prefer scripting, open the default script *matRad.m* from the main matRad folder

.. code-block:: matlab

    edit matRad.m

You can use it to learn something about the code structure of a typical matRad workflow and execute it section by section.

If you want to run the full script for an exemplary photon plan, you can type

.. code-block:: matlab

    matRad

in your command window.

.. rubric:: Option 3: Using the examples
    :heading-level: 2

This is by far the most time consuming, but also the most educational approach to matRad, as it includes working through tutorial-style example scripts for different scenarios, e.g. treatment planning with different radiation sources or including biological models. Therefore, if you want to have in-depth training on the whole bandwidth of applications of our software, working through the examples is what we recommend.

In order to get started with the examples, you simply need to navigate to the folder *examples* when in the main matRad directory, open any example of your choice, execute it section by section and move on to the next example afterwards. As the code includes detailed instructional comments, the examples serve as a step-by-step guide to different matRad use cases and are easy to follow and understand.

**See also:**

.. toctree::
   :maxdepth: 1

   guiintro
   matradscript

..

