.. _basedata_particles:

.. toctree::
   :maxdepth: 2

========================
Particle Base Data File
========================

This page describes the format and the variables as part of the particle base data files `protons_Generic.mat <https://github.com/e0404/matRad/blob/master/basedata/protons_Generic.mat>`_ and `carbon_Generic.mat <https://github.com/e0404/matRad/blob/master/basedata/carbon_Generic.mat>`_. 
The base data is stored using MATLAB's structure format. 
The first sub level of the structure 'machine' contains two fields named 'meta' and 'data' which are explained separately next.

machine.meta
------------
Stores relevant isolated meta information describing the actual base data file that does not fit into the machine.data structure.

machine.meta.machine
^^^^^^^^^^^^^^^^^^^^
Actual machine name as string to identify the base data set during runtime.

machine.meta.radiationMode
^^^^^^^^^^^^^^^^^^^^^^^^^^
The radiationMode stores the radiation modality the base data is describing. Logically, `protons_Generic.mat <https://github.com/e0404/matRad/blob/master/basedata/protons_Generic.mat>`_ models protons in water and the corresponding field machine.meta.radiationMode is set to 'protons'.

machine.meta.SAD
^^^^^^^^^^^^^^^^
This subfield holds the geometrical source to axis distance in millimeter. In case of the generic base data set, we use a value of 10000 [mm].

machine.meta.BAMStoIsoDist
^^^^^^^^^^^^^^^^^^^^^^^^^^
This subfield depicts the geometrical distance from the beam application monitoring system/beam nozzle to the isocenter. For the generic base data set we use a value of 2000 [mm].

machine.meta.LUT_bxWidthminFWHM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This field is a two-dimensional array and holds a look up table which describes a mapping of the lateral spot spacing (bixel width) to the focus index (initial beam width). In the generic proton base data set, we only store one focus index (initial beam width) for which reasons the corresponding LUT_bxWidthminFWHM is set to [1,Inf;5 5]. This means that if the lateral spot spacing is set during treatment planning to a value between 1 and Inf, then a focus index with at least a full width half maximum (FWHM) of 5 mm is consistently used. As we only store one focus index in machine.data.initFocus for `protons_Generic.mat <https://github.com/e0404/matRad/blob/master/basedata/protons_Generic.mat>`_ with a FWHM at isocenter for the highest beam energy of ~5.4 [mm], we basically make sure to always use the first (and only) focus index. In principle, the idea behind machine.meta.LUT_bxWidthminFWHM is to allow for different focus indices when using different lateral spot spacings (bixel widths). Users can model the following behavior: E.g. when setting the lateral spot spacing to 5 [mm] then use a foci index with a FWHM greater than 8 [mm]. Or when setting the lateral spot spacing to 3 [mm] then use a focus index with FWHM greater than 6 [mm].  

machine.data
------------
This subfield contains an array of structures holding the corresponding depth-dependent physical and (biological) beam properties for each initial beam energy.

machine.data.energy
^^^^^^^^^^^^^^^^^^^
Initial beam energy in MeV/u.

machine.data.depths
^^^^^^^^^^^^^^^^^^^
Depth values stored on an irregular grid in [mm] (higher resolution around the peak). All depth-dependent data like (Z, sigma) are stored exactly on these depth positions.

machine.data.peakPos
^^^^^^^^^^^^^^^^^^^
Peak position of the pencil beam in [mm].

machine.data.offset
^^^^^^^^^^^^^^^^^^^
This field allows to consider a pencil beam offset caused by passive beam line elements, that have not been modeled during the creation of the base data set. 

machine.data.Z
^^^^^^^^^^^^^^
This field holds the integrated depth dose profiles of the corresponding radiation modality. Regarding units, we refer to the base data section in the section :ref:`Dose influence matrix calculation <dose_calc>`.

machine.data.sigma 
^^^^^^^^^^^^^^^^^^^
(for a single Gaussian lateral beam model) 
Lateral beam broadening of the particle pencil beam in water. Hint: The first sigma value for depth 0 [mm] should be set to 0 [mm] because the initial beam width at the patient surface is modeled via machine.data.initFocus and covered later in this wiki page.

machine.data.sigma1 
^^^^^^^^^^^^^^^^^^^
(for double Gaussian lateral beam model) 
Lateral beam broadening of the particle pencil beam in water of the narrow Gaussian component. Again, the first sigma1 value for depth 0 [mm] should be set to 0 [mm].

machine.data.sigma2 
^^^^^^^^^^^^^^^^^^^
(for double Gaussian lateral beam model) 
Lateral beam broadening of the particle pencil beam in water of the broad Gaussian component.

machine.data.w
^^^^^^^^^^^^^^ 
(for double Gaussian lateral beam model) 
Relative weight between the narrow (sigma1) and the broad (sigma) Gaussian component.

machine.data.initFocus
^^^^^^^^^^^^^^^^^^^^^^
Let numFoci be the number of available focus indices and machine.data.initFocus hold three subfields named 'dist', 'sigma' and 'SisFWHMAtIso' of the following dimensions numFoci x N, numFoci x N and numFoci x 1 whereas N indicates the number of values used in the look up table. SisFWHMAtIso describes for each focus index the initial FWHM at isocenter. In contrast, 'dist' and 'sigma' depict a look up table to model the particle beam spread in air. In the generic proton and carbon base data set we do not model beam widening in air from the beam nozzle to the patient surface (although the code is capable of). Therefore, machine.data.initFocus(1).sigma is a constant value over a distance from 0 to 20000 [mm].
