.. _visualization:

=============
Visualization
=============

Graphical User Interface
------------------------

At any point, it is possible to start the graphical user interface by entering 

.. code-block:: matlab
    
    matRadGUI

in the command window. The GUI will automatically check your base workspace and recognizes your current planning step. Data will be displayed, according to the detected planning step. Please note that a dose cube will be displayed not until you have finished an optimization.

It is possible to adjust the plot to your needs in the lower left area of the GUI. For instance, it is possible to switch to a longitudinal profile plot of the central beam axis.

Plotting tools
--------------

There's a number of tools to visualize results from the command line or per script.

* The top-level function :func:`matRad_planAnalysis` provides a convenient way to visualize the dose distribution and DVHs of a plan. It can be used in the command line or in scripts, and it automatically detects the current planning step and displays the relevant data.
* :class:`matRad_ViewingWidget` only opens the Viewer part of the GUI for quick visualization.
* :func:`matRad_plotSlice` can be used to plot a single slice of the dose distribution with underlying CT and many configuration options.
* :func:`matRad_showDVH` can be used to plot DVHs of a plan.


