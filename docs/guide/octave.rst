.. _octave:

=============================
How to run matRad with Octave
=============================

We try to keep matRad's core compatible with `GNU Octave <https://www.gnu.org/software/octave/>`_ to allow users to run matRad without a MATLAB license. Even the graphical user interface (GUI) can, with limitations, be used. To maintain compatibility, nearly all tests are currently run on GNU Octave 6 in addition to MATLAB.

One difficulty with running matRad in GNU Octave is the optimization module relying on `Ipopt <https://github.com/coin-or/Ipopt>`_ with a MEX interface.

While we provide some GNU Octave mex files for a few Windows versions, the IPOPT mex files usually needs to  be compiled from source with Octave. The steps to compile the Ipopt interface with Octave are presented below.

Installing GNU Octave
---------------------

It is recommended to use matRad with the latest stable release of `GNU Octave <https://www.gnu.org/software/octave/>`_. matRad has been mainly tested with GNU Octave versions 6 in Linux, but also newer versions seem to work. When installing Octave from the package manager of a Linux distribution, it is necessary to also install the development files to provide the command ``mkoctfile`` required to build the interface to Ipopt.

By default, Octave is distributed with 32-bit indexing on Linux, while for Windows also a 64bit indexing installer is available. GNU Octave can also be compiled from source with 64-bit indexing. However, building Octave with 64-bit indexing is rather laborious requiring the compilation of several external libraries to enable 64-bit indexing. `Click here <https://www.gnu.org/software/octave/doc/v4.2.0/Compiling-Octave-with-64_002dbit-Indexing.html>`_ for further details on how to compile GNU Octave with 64-bit indexing. Notice that the 64-bit indexing option of Octave is still experimental and that special care should be taken when linking the libraries to avoid segmentation faults.

While matRad can be used with Octave compiled either with 32-bit or 64-bit indexing, some clinically relevant treatment planning scenarios may require an influence matrix, which exceeds the maximum array size available with 32-bit indexing.

IPOPT
-----

Since version 3.1.0, the `IPOPT folder <https://github.com/e0404/matRad/tree/7df8436bbcb881fb09a3aafd94f2ef39cd0dfb66/thirdParty/IPOPT>`_ in matRad contains precompiled mex files for certain Octave versions in Windows. While in MATLAB, mex files have an extension naming scheme depending on the operating system and are often compatible across multiple MATLAB versions, we have oberserved that this is not true for GNU Octave, where a mex file just always uses "mex" as extension. Thus matRad uses a dedicated naming scheme to store Octave mex files as "mexoct" + version + os + arch. For Octave 6.4.0 on Windows, for example, the filename would be ipopt.mexoct640w64. When running IPOPT from Octave, matRad will copy this file to ipopt.mex to call it.

The MATLAB interface of Ipopt distributed with matRad is linked against a MATLAB-specific library for the linear solver MA57. An equivalent interface of Ipopt for GNU Octave should be compiled from source. It requires the installation of the `Ipopt <https://github.com/coin-or/Ipopt>`_ library along with the corresponding header files. Since the linear solvers are not distributed along with Ipopt, they should be obtained from third-party software. Different linear solvers libraries can be used including `HSL_MA57 <http://www.hsl.rl.ac.uk/catalogue/ma57.html>`_, `MUMPS <http://mumps.enseeiht.fr/>`_ and others. In addition to compiling Ipopt and the linear solvers from source, Ipopt linked to MUMPS is typically available from the package manager of the Linux distributions. Ipopt version 3.12.8 linked with the linear solvers HSL_MA57 and MUMPS has been successfully applied in matRad with GNU Octave.

Compiling Octave interface of Ipopt
-----------------------------------

The source code of Ipopt provides the source files to compile the original MATLAB interface of Ipopt distributed with matRad. In order to compile the equivalent interface for GNU Octave, it is recommended to use the `interface rewritten by Enrico Bertolazzi <https://github.com/ebertolazzi/mexIPOPT>`_. The rewritten interface requires GNU Octave version 4.2 or above. The Octave interface of Ipopt can be compiled using the ``mkoctfile`` command with the option to produce MEX files. 

Since version 3.1.0, matRad provides scripts to help with compilation of the Ipopt interface for GNU Octave. The scripts target Windows platforms and can be executed from within the MinGW environment installed with GNU Octave by using the ``cmdshell.bat`` in Octave's root folder. These scripts are located in the `thirdParty/IPOPT` folder and can also be easily adapted by linux users.