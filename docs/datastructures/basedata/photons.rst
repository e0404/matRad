.. _basedata_photons:

======================
Photon Base Data File
======================

This page describes the format and the variables as part of the particle base data files `protons_Generic.mat <https://github.com/e0404/matRad/blob/master/basedata/protons_Generic.mat>`_ and `carbon_Generic.mat <https://github.com/e0404/matRad/blob/master/basedata/carbon_Generic.mat>`_.
The base data is stored using MATLAB's structure format.
The first sub level of the structure 'machine' contains two fields named 'meta' and 'data' which are explained separately next.

machine.version
---------------
Specifies the version of the machine file. If not present, version 1 is assumed. Currently version 1 and 2 are supported. 
Version 2 is the current standard used in matRad. It's main difference to version 1 is that kernel normalization does not imply a convolution grid of 0.5 mm.

machine.meta
------------
Stores relevant isolated meta information describing the actual base data file that does not fit into the machine.data structure.

machine.meta.name
^^^^^^^^^^^^^^^^^
Actual machine name as string to identify the base data set during runtime.

machine.meta.radiationMode
^^^^^^^^^^^^^^^^^^^^^^^^^^
The radiationMode stores the radiation modality the base data is describing. Logically, `protons_Generic.mat <https://github.com/e0404/matRad/blob/master/basedata/protons_Generic.mat>`_ models protons in water and the corresponding field machine.meta.radiationMode is set to 'protons'.

machine.meta.SAD
^^^^^^^^^^^^^^^^
This subfield holds the geometrical source to axis distance in millimeter. In case of the generic base data set, we use a value of 1000 [mm].

machine.meta.SCD
^^^^^^^^^^^^^^^^
This subfield holds the geometrical source to collimator distance in millimeter. In case of the generic base data set, we use a value of 500 [mm].
This value is mainly used to convert between penumbra and source width.

machine.data
------------
This subfield contains the detailed dosimetric machine information. It is currently only storing decomposed pencil-beam kernels according to `Bortfeld et al. (1993) Medical Physics <https://doi.org/10.1118/1.597070>`_.

machine.data.energy
^^^^^^^^^^^^^^^^^^^
Maximum photon energy in MeV. Corresponds to the acceleration potential in MV. 

machine.data.betas
^^^^^^^^^^^^^^^^^^
Beta factors for the individual kernels. Row array with length of number of kernels.

machine.data.m
^^^^^^^^^^^^^^
Attenuation factor (obtained from kernel fitting).

machine.data.primaryFluence
^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
Radial fluence distribution of a photon beam.

machine.data.kernelPos
^^^^^^^^^^^^^^^^^^^^^^
Lateral Kernel evaluation positions in [mm]. Alle kernels in machine.data.kernel are tabulated at these positions.

machine.data.penumbraFWHMatIso
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Penumbra full width half maximum (FWHM) at isocenter in [mm]. 
This value is used to convert between penumbra and source width and in kernel convolution with the field function.

machine.data.kernel
^^^^^^^^^^^^^^^^^^^
Structure Array with all kernels by SSD. 
Each kernel struct contains the field ``SSD`` to tell the corresponding SSD for which the kernel components were fitted.
Further fields enumerate the kernels, i.e.: ``kernel1``, ``kernel2``, ``kernel3``, etc. 
Each kernel is tabulated at the positions given in ``machine.data.kernelPos`` and has been fitted to the measured data.
Number of kernels must match the length of ``machine.data.betas``.
