.. _pln:

======================
The pln Data Structure
======================

The ``pln`` struct holds the meta information about the radiation treatment plan.

Top-level properties stored in the ``pln`` struct
-------------------------------------------------

Information relevant for all parts of the treatment planning workflow is stored on the top-level of the struct.
These include the radiation modality, the machine information, the number of fractions, and scenario models as well as biological models.

**pln.radiationMode**

    Specifies the radiation modality. Can either be *photons*, *protons*, *helium* or *carbon*.

**pln.machine**

    In order to load the appropriate base data, one can introduce the machine in which the treatment plan has been incorporated and the code will look for the file under ``{pln.radiationMode}_{pln.machine}.mat``. For example, setting ``pln.machine`` to 'Generic' for a photon treatment plan will load the already available ``photons_Generic.mat`` file.

**pln.numOfFractions**

    Specifies the number of fractions. Note that this parameter only needs to be set for biological treatment planning for carbon ions, where the optimization process is based on the fraction dose and not based on the overall dose.

**pln.multScen**
    Specifies a scenario model for the treatment plan, see :mod:`matRad.scenarios`.

**pln.bioModel**
    Specifies a biological model for the treatment plan, see :mod:`matRad.bioModels`.

Workflow Step Configuration properties
--------------------------------------

Workflow step configuration properties can be stored in the ``pln`` struct.
The Syntax for accessing these properties is ``pln.prop{StepName}.{PropertyName}``. 
This results in a nested structure, where the first level is the step name and the second level is the property name.

Current possible names are  ``propStf`` (steering information / geometry), ``propDoseCalc`` (dose calculation), ``propOpt`` (optimization), and ``propSeq`` (sequencing).


.. admonition:: Using the ``pln.prop{StepName}.{PropertyName}`` vs direct class instantiation.
   :class: note

   matRad distinguishes between a top-level API and low-level programming using the classes defining workflow steps.
   Using top-level functions like :func:`matRad_calcDoseInfluence` from the root :mod:`matRad` folder will take pln as an argument, instantiate the appropriate object (e.g. :class:`matRad_ParticleHongPencilBeamEngine`), and try to configure its properties from the ``pln.prop{StepName}.{PropertyName}`` structure.
   Alternatively, these classes can be direclty used and properties can be explicitly set.

Here's an overview of the property mapping:

.. include:: ../includes/planapi.rst

pln.propStf
^^^^^^^^^^^

**pln.propStf.gantryAngles**

    Specifies the gantry angles as MATLAB vector according to the `matRad coordinate system <The-matRad-coordinate-system>`_.

**pln.propStf.couchAngles**

    Specifies the couch angles as MATLAB vector according to the `matRad coordinate system <The-matRad-coordinate-system>`_.

**pln.propStf.bixelWidth**

    Specifies the width (and height) of quadratic photon bixels (i.e. discrete fluence elements). For particles, this parameter specifies the lateral spot distance.

**pln.propStf.numOfBeams**

    Specifies the number of beam directions. During the matRad script, this parameter is automatically determined.

**pln.propStf.isocenter**

    Specifies the isocenter of the treatment plan in voxel coordinates within the ct.cube. By default, the isocenter is calculated as the center of gravity of all voxels belonging to structures that have been modeled as target volume in the `cst cell <The-cst-cell>`_.

pln.propOpt
-----------

**pln.propOpt.bioOptimization**

    Specifies the type of biological optimization. *none* corresponds to a conventional optimization based on the physical dose. *effect* corresponds to an effect based optimization according to `Wilkens & Oelfke <http://iopscience.iop.org/0031-9155/51/12/009>`_. *RBExD* corresponds to an optimization of the RBE weighted dose according to `Krämer & Scholz <http://iopscience.iop.org/0031-9155/51/8/001>`_.

**pln.propOpt.runDAO**

    Setting this value to ``true`` will enable Direct Aperture Optimization run allowing us to directly optimize aperture shapes and weights.

**pln.propOpt.runSequencing**

    Setting this value to ``true`` will enable sequencing algorithms run.

Additional adjustable properties
-------------------------------

The following properties of the pln struct can additionally be adjusted. If they are not explicitly set, default values are used. The default values are handled by the `MatRad_Config class <https://github.com/e0404/matRad/blob/master/MatRad_Config.m>`_.

**pln.propStf.longitudinalSpotSpacing**

    Specifies the longitudinal spot spacing. Default: *3* mm.

**pln.propStf.addMargin**

    If this property is set to *true*, the target is expanded for beamlet finding. Default: *true*.

**pln.propDoseCalc.doseGrid.resolution**

    Specifies the resolution for the dose calculation. Default: x direction: *3* mm, y direction: *3* mm, z direction: *3* mm.

**pln.propDoseCalc.defaultLateralCutOff**

    Specifies the lateral cutoff. Default: *0.995* rel.

**pln.propDoseCalc.defaultGeometricCutOff**

    Specifies the geometric cutoff. Default: *50* mm.

**pln.propDoseCalc.ssdDensityThreshold**

    Specifies the ssd density threshold. Default: *0.05* rel.

**pln.propOpt.defaultMaxIter**

    Specifies the number of maximum iterations. Default: *500*.

**pln.propMC.ompMC_defaultHistories**

    Specifies the number of particles simulated per pencil beam photon Monte Carlo calculation ompMC. Default: *1e6*.

**pln.propMC.ompMC_outputVariance**

    If it is set to *true*, variance scoring for MCsquare is not supported for the photon Monte Carlo calculation ompMC. Default: *false*.

**pln.propMC.MCsquare_defaultHistories**

    Specifies the number of particles simulated per pencil beam for the proton Monte Carlo calculation MCsquare. Default: *1e6*.

**pln.propMC.direct_defaultHistories**

    Specifies the number of particles simulated per pencil beam for the Monte Carlo calculation when bypassing the dij calculation. Default: *2e4*.

**pln.disableGUI**

    If this value is set to *true*, matRadGUI is disabled. Default: *false*.