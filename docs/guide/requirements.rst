.. _requirements:

===========================
Minimum System Requirements
===========================

Software Environment
--------------------

MATLAB
^^^^^^

We develop and test matRad mainly with MATLAB version **2022b and later**. Our automatic testing framework via GitHub Actions currently uses R2022b as well as the latest MATLAB version. If you run into a problem with these or any other version, please let us know, but do not expect us to try to stay compatible with all other versions. The main reason for limiting us to these MATLAB versions are incompatibilities in the mex interface to IPOPT (especially on Windows). Other incompatibilities often reveal themselves as missing functions (like ``isstring``, for example).

MATLAB Toolboxes
^^^^^^^^^^^^^^^^

To run all dose calculation and optimization functionalities you will only need a basic MATLAB installation. For DICOM import and export, you also need the Image Processing Toolbox. Certain additional functionalities are available with the Optimization Toolbox (``fmincon``), the Parallel Computing Toolbox, and the Statistics & Machine Learning Toolbox.

Octave
^^^^^^

Note that compatibility with Octave is not our primary goal, but it is also part of the automatic testing framework on GitHub Actions for Octave 6.4.

Operating System
----------------

Since we work in the programming environment MATLAB, operating system incompatibilities are not that common. They may arise, in particular, when using mex interfaces. Our precompiled mex interfaces should work on Windows 10 & 11, Ubuntu 18.04 (and later), and MacOS High Sierra. Please let us know if you run into issues, but the first step should always be trying to compile the mex interfaces yourself.

Especially on Mac, there might be substantial issues due to their annoying quarantine system, which prevents the execution of downloaded files. If you run into this problem, you need to remove the respective quarantine flags, especially from mex files you'd like to use.

Hardware Requirements
---------------------

There are no hard minimum requirements to do dose calculation and optimization with matRad. We do treatment planning tutorials also with systems with 2GB RAM, but that means that the cases you are looking at are somewhat small (low spatial resolution, few beams, rather no particles). If you want to do treatment planning at realistic resolutions, we recommend 16GB RAM or more.

If you run into memory problems, you have basically three options:

* Buy more RAM ;)
* Import your data at low spatial resolution, which is possible during DICOM import. Remember that reducing the resolution by a factor of 2 will reduce memory consumption by a factor of 8! Alternatively or additionally, reduce the resolution of the dose calculation grid using the option ``pln.propDoseCalc.doseGrid.resolution``.
* Reduce the number of pencil-beams by choosing larger ``bixelWidth`` or, for particles only, the longitudinal spot spacing.
* Restrict the dose calculation area by specifying tighter lateral cut-off values in `matRad_calcParticleDose (line 211) <https://github.com/e0404/matRad/blob/master/matRad_calcParticleDose.m#L211>`_ and `matRad_calcPhotonDose (line 61) <https://github.com/e0404/matRad/blob/master/matRad_calcPhotonDose.m#L61>`_, respectively. While this induces inaccuracy in the planning process, this might be a viable option for educational purposes.