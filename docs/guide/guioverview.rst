.. _guioverview:

##############################################
Overview of the graphical user interface (GUI)
##############################################

The GUI is used for the treatment plan visualization and the adjustment of the plan and optimization parameters. It is possible to execute the matRad script using only the GUI (see :ref:`run_gui`).

Overview of the GUI
===================

The matRad GUI consists of 6 main sections:

.. list-table::
    :header-rows: 1

    * - Section
      - Content
    * - :ref:`Workflow <workflow>`
      - Includes the main steps that need to be executed for treatment planning.
    * - :ref:`Plan <planParameters>`
      - Here the plan parameters like beam direction and radiation mode can be selected.
    * - Command Window
      - Here the output of the Command Window is displayed inside the GUI.
    * - :ref:`Objectives & constraints <optParameters>`
      - Using the optimization parameters the constraints of the VOIs can be adjusted.
    * - :ref:`Visualization <visualization>`
      - Here the visualization can be adjusted to show different planes/slices of plot types.
    * - Viewing
      - The information specified in the Visualization Parameters is displayed in this section.

.. image:: /images/GUI-Guide_optimizedGUIScreenshot.png
    :width: 650px

.. _workflow:

Workflow
========

In the Workflow section, the patient data is initially loaded. You can also start the :ref:`dicom import <dicomimport>` from here. After the adjustment of all parameters, the dose calculation and the fluence optimization can be started from here:

.. image:: /images/GUI-Guide_workflowGUIScreenshot.png

.. _planParameters:

Adjustment of the plan parameters
=================================

In this section, the plan parameters are adjusted before calculating the dose-influence-matrix.

.. list-table::
    :header-rows: 1

    * - Parameter
      - Description
    * - bixel width
      - *Photons:* width of a photon bixel. *Particles:* lateral spot distance. Default value: 5 mm
    * - Gantry angle
      - Here the set of desired gantry angles (in degree) can be specified. For the separation of the values you can either use ',' or a space ' '.
    * - Couch angle
      - Here the set of desired couch angles (in degree) can be specified. For the separation of the values you can either use ',' or a space ' '. Make sure that you always have the same amount of couch and gantry angles.
    * - Radiation mode
      - You can choose between photons, protons and carbon ions.
    * - Machine
      - Actual machine name as string to identify the base data set during runtime.
    * - IsoCenter
      - Use this to set the isocenter (in mm).
    * - # Fractions
      - Here the desired number of fractions can be specified.
    * - Biological optimization
      - For carbon ions, you can apply a biological optimization. You can choose between an optimization of the biological effect (``effect``) or the RBE-weighted dose (``RBExD``).
    * - Run Sequencing
      - Check this if you want to run a MLC sequencing. The number of stratification levels can be adjusted.
    * - Run Direct Aperture Optimization
      - Check this if you want to run an additional direct aperture optimization.

.. image:: /images/GUI-Guide_planParametersGUIScreenshot.png

.. _optParameters:

Adjustment of the optimization parameters
=========================================

The optimization parameters regarding the volumes of interest (VOIs) are stored in the variable ``cst``. For more detailed information about the parameters stored in the cell, please refer to the :ref:`documentation of the cst-cell <cst>`. Using the GUI, you can adjust the settings. To add or delete volumes, you can use the "+" and "-" buttons.

.. list-table::
    :header-rows: 1

    * - Field
      - Description
    * - VOI name
      - Via a drop-down menu, you can select a VOI by clicking its name.
    * - VOI type
      - You can specify whether the VOI is an organ at risk (OAR) or a target volume.
    * - OP
      - *Overlap*. This value defines how overlapping structures are handled during optimization. Consider two structures A and B with priorities :math:`p_A` and :math:`p_B`. If A and B both include voxel *i*, voxel *i* will be treated to belong only to structure A if :math:`p_A < p_B`. If :math:`p_A = p_B` the voxel will be considered for both structures. An extension to more than two structures is trivial.
    * - Function
      - *Objective Function*. This field allows you to specify how the VOI will be considered during the optimization. You can choose between *Squared Underdosing*, *Squared Overdosing*, *Squared Deviation*, *Mean Dose*, *EUD*, *Max DVH*, *Min DVH*, *DVH constraint*, *Min/Max dose constraint* or *mean dose constraint*. You can find more detailed information about this in the section 'Dose objectives' of the page: :ref:`The cst cell <cst>`.
    * - p
      - *Penalty*. For the objective function value, a weighted sum is calculated. The penalty value corresponds to the weighting factor for this VOI with respect to the defined constraint (e.g. overdosing). By adjusting this value, you can stress the importance of these constraints with respect to each other.
    * - Parameters
      - For *Squared Underdosing*, *Squared Overdosing* and *Squared Deviation* this value corresponds to the threshold dose above/below which the penalty will apply. For the *Mean Dose* option, this value is not needed, as the mean dose within this VOI will be minimized. For the *EUD* method, the parameter corresponds to the exponent.

.. _visualization:

Visualizing treatment plans
==========================

After the optimization, the treatment plan can be visualized within the GUI. Using the visualization parameters, you can change the view. The radio buttons can be used to turn on or off, among others, the plotting of contours, dose (isolines), and isoline labels.

.. image:: /images/doseVisParameter.png

Display options
---------------

.. list-table::
    :header-rows: 1

    * - Type of plot
      - Display option
      - Plane
      - Resulting image
    * - **intensity**
      - **Dose**
      - axial
      - .. image:: /images/doseVisAxialIntensity.png
    * -
      -
      - sagittal
      - .. image:: /images/doseVisSagitalIntensity.png
    * -
      -
      - coronal
      - .. image:: /images/doseVisCoronalIntensity.png
    * -
      - **effect**
      - axial
      - .. image:: /images/doseVisAxialEffect.png
    * -
      - **RBEWeightedDose**
      - axial
      - .. image:: /images/doseVisAxialRBExD.png
    * -
      - **RBE**
      - axial
      - .. image:: /images/doseVisAxialRBE.png
    * -
      - **alpha**
      - axial
      - .. image:: /images/doseVisAxialAlpha.png
    * -
      - **beta**
      - axial
      - .. image:: /images/doseVisAxialBeta.png
    * -
      - **RBETruncated10Perc**
      - axial
      - .. image:: /images/doseVisAxialRBEtruncated.png
    * - **profile**
      -
      - **lateral**
      - .. image:: /images/doseVisLateralProfile.png
    * -
      -
      - **longitudinal**
      - .. image:: /images/doseVisLongitudinalProfile.png

DVH
---

To draw a DVH of the current treatment plan and display some quality indicators you can click the *Show DVH/QI* button:

.. image:: /images/DVHVisScreenshot.png