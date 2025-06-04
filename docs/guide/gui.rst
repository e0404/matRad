.. _run_gui:

##########################################
Running the graphical user interface (GUI)
##########################################

To execute the matRad GUI using MATLAB you need to:

.. contents::
   :local:
   :depth: 2

If you prefer to use the :file:`matRad.m` script to execute matRad, check out the :ref:`matRad script <run_script>`.  
For more detailed information about the different features of the GUI you can take a look at :ref:`matRad GUI Overview <guioverview>`.

Step 1: Open matRad folder in MATLAB
------------------------------------

To use matRad you need to open the matRad folder in MATLAB.
Open MATLAB and navigate to the location of the files, if you have cloned the repository it is most likely located in your local Github folder.

Inside the :ref:`root folder <root>` of the matRad repository there are only a few files, of which the most important ones are::

* :func:`matRad_rc.m <matRad_rc>` - the matRad configuration script setting path and environment
* :scpt:`matRad.m <matRad>` - the main introductory script to run matRad workflow
* :func:`matRadGUI.m <matRadGUI>` - the main script to run the matRad GUI

The function to run the matRad GUI is called :func:`matRadGUI`, which instantiates the :class:`matRad_MainGUI`.

Step 2: Start the matRad GUI
----------------------------

To start the GUI select :file:`matRadGUI.m` from your current folder and run it (right-click → run or F9) or simply type ``matRadGUI`` in your command window.  
Now the empty GUI should be opened:

.. image:: /images/GUI-Guide_emptyGUIScreenshot.png
    :width: 650px

If the GUI is not empty, then there is a patient already loaded in your workspace. To get an empty GUI you can clear your workspace and restart the GUI. However, this is not necessary as you can simply load a new patient.

Step 3: Execute treatment planning
----------------------------------

**Load patient data**

First, you need to load the patient data. Therefore, the matRad release contains the :doc:`The-CORT-dataset`. To import additional patient data have a look at the :doc:`The-dicom-import`.  
To load a patient click the **Load \*.mat data** button in the **Workflow** section.  
A window should open. In the folder ``phantoms``, you can find different patient files.

.. image:: /images/GUI-Guide_loadDataGUIScreenshot.png

Here you can select which patient file (``*.mat``) you want to load. Upon opening the ``*.mat`` file the patient data is loaded into the GUI:  
On the right side of the GUI you should see the patient-CT with the defined VOIs. On the left side, the optimization parameter table should now be filled.

.. image:: /images/GUI-Guide_loadedGUIScreenshot.png
    :width: 650px

**Set plan parameters**

Now you can start to adjust the plan parameters:

+-----------------------------------------------------------------------------------------------+
| The bixel width, as well as the isocenter, can be adjusted but should already be set to       |
| reasonable values.                                                                            |
+-----------------------------------------------------------------------------------------------+
| To set the beam directions you have to select the according gantry and couch angles. Every    |
| gantry angle defines a beam and needs a couch angle.                                          |
+-----------------------------------------------------------------------------------------------+
| For the radiation mode, you can choose photons, protons or carbon.                            |
+-----------------------------------------------------------------------------------------------+
| If you set carbon as radiation mode, you can activate the biological optimization. You can    |
| choose between an effect based optimization (*effect*) or the optimization of the RBE-weighted|
| dose (*RBExD*).                                                                               |
+-----------------------------------------------------------------------------------------------+
| For the radiation mode "photons", you have the option to run a MLC sequencing, where you can  |
| set the number of stratification levels and additionally you can run a direct aperture        |
| optimization.                                                                                 |
+-----------------------------------------------------------------------------------------------+

.. image:: /images/GUI-Guide_planParametersGUIScreenshot.png

**Set optimization parameters**

The optimization parameters are used to influence the outcome of the fluence optimization. Here you can set the parameters of the VOIs (e.g. min/max dose, penalty, overlap priority, etc.). For more information, take a look at the :doc:`The-cst-cell`. Using the '**+**' and '**-**' buttons you can add and remove VOIs.

The column ``p`` (*penalty*) determines the relative weighting of the objective within the overall weighted sum objective function. The column ``Parameters`` lets you specify additional parameters for given objectives. For squared over- and underdosage as well as squared deviation, this simply corresponds to the reference dose level, for EUD it is the exponent. A mean dose objective does not require an additional parameter.

The column ``OP`` (*overlap priority*) is very important as it determines the assignment of voxels to structures. Consider a voxel that belongs to two structures, e.g. to the rectum and to the prostate. For the optimizer it is necessary to distinguish to which structure the voxel should belong to during optimization. If you assign priority 1 to the prostate and priority 2 to the rectum in our example, every voxel within the overlap of the two structures will be considered as prostate during optimization. If you assign priority 2 to the prostate and priority 1 to the rectum in our example, every voxel within the overlap of the two structures will be considered as rectum during optimization. Assigning the same priority to overlapping structures will result in the overlapping volume being considered for all structures. Be aware that the skin contour usually encompasses the entire patient volume. If you want to make sure that target voxels are not also considered skin voxels you need to assign a lower priority (i.e. a higher number) to the skin volume. Please check with the provided example patient data to understand this in full detail.

*Note: Changing the VOI Type from OAR to target will lead to additional beamlets or spots that need to be considered for the dose-influence-matrix calculation. As a result, these changes have to be done before the Dij-calculation.*

.. image:: /images/GUI-Guide_optimizationParametersGUIScreenshot.png

**Calculate Dose influence matrix**

To start the calculation of the dose-influence-matrix you simply need to click the **Calc. Dose Influence** button in the workflow:

.. image:: /images/GUI-Guide_workflowGUIScreenshot.png

You should see a window pop up, showing a progress bar of the calculation:

.. image:: /images/GUI-Guide_dijProgressBarScreenshot.png

In addition, the progress is displayed in the Command Window:

.. image:: /images/GUI-Guide_dijOutputScreenshot.png

**Execute fluence optimization**

Once the dose calculation is completed, you can start the fluence optimization by clicking the **Optimize** button in the workflow section. The iterations of the optimization are displayed in the Command Window:

.. image:: /images/GUI-Guide_fluenceOptOutputScreenshot.png

To adjust the convergence criteria you can specify the *maximum number of iterations* and the *convergence* precision in the *Optimization Parameter* section. Default values are: 1000 iterations and a precision of :math:`10^{-3}`:  
(Precision ≡ |(FuncValue_old − FuncValue_new) / FuncValue_old|)

.. image:: /images/GUI-Guide_optimizationParameters2.png

Step 4: Visualize resulting treatment plan
------------------------------------------

Once the fluence optimization has converged the resulting dose distribution will be displayed in the GUI. Here you can adjust the visualization parameters to display different slices/planes, use different plot types, etc.

.. image:: /images/GUI-Guide_optimizedGUIScreenshot.png
    :width: 650px

To calculate a DVH of all VOIs and to see the quality indicators (which contain the mean/max/min dose for each VOI) you can use the **Show DVH/QI** button in the *Visualization* section.

.. image:: /images/DVHVisScreenshot.png


Step 5: Import additional patient data
--------------------------------------

matRad supports the import of patient data stored in the DICOM format. A set of functions designed for this purpose can be found in the subfolder :file:`dicom <https://github.com/e0404/matRad/tree/master/dicom>`. For more information about the usage of the import functions please check out :doc:`The-dicom-import`.

