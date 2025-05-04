.. _faq:

.. toctree::
   :maxdepth: 2

==========================
Frequently Asked Questions
==========================

While this site attempts to cover a set of frequently asked questions, it is not exhaustive.
We suggest to visit our `GitHub Discussion Forum <https://github.com/e0404/matRad/discussions>`_ for Questions and Answers from the Community.

.. admonition:: How can I model custom particle machines in matRad?
    :class: dropdown

    You need to provide your own base data file according to our examples (*proton_Generic.mat*, *carbon_Generic.mat*). More information about the format is summarized on :ref:`this page <basedata_particles>`.


.. admonition:: How can I use a custom HLUT table in matRad?
    :class: dropdown

    matRad has the following procedure for reading HU lookup tables via `matRad_loadHLUT.m <https://github.com/e0404/matRad/blob/c22da7d2a593b7ae56277ce124074c55fe4d45dd/dicom/matRad_loadHLUT.m>`_ during the dose calculation:
    *  matRad first checks if there is any `.hlut` file provided by the user in the `hlutlibrary <https://github.com/e0404/matRad/tree/master/dicom/hlutLibrary>`_ directory, with the following naming convention: MANUFACTURER-MODEL-ConvolutionKernel_CONVOLUTIONKERNEL_RADIATIONMODALITY.hlut (e.g. Philips-AcQSimCT-ConvolutionKernel-000000_protons.hlut). 
    *  If there is no file with this name available, matRad would use its own default lookup table with the following naming convention: matRad_default_RADIATIONMODALITY.hlut (e.g. `matRad_default.hlut <https://github.com/e0404/matRad/blob/c22da7d2a593b7ae56277ce124074c55fe4d45dd/dicom/hlutLibrary/matRad_default.hlut>`_). This file comprises two columns, first is the HU units and second is the corresponding electron density/stopping power (comments indicated by starting a line with "#" will be omitted).

    Generating your own custom HLUT table can therefore be done in two ways, either making a custom file with the mentioned naming convention or changing the default tables provided by matRad (not recommended).

.. admonition:: Why is the matRad GUI slow (macOS)?
    :class: dropdown

    This issue may be caused by conflicts between MATLAB and window-snapping apps on macOS (Rectangle, BetterTouchTool, Magnet etc.). 
    May be remedied by quitting the app or disabling the "Window Snapping " feature then restarting MATLAB.
    For further information refer to : 
    * `matRad Issue #550 <https://github.com/e0404/matRad/issues/550>_`
    * `MATLAB Answers thread <https://de.mathworks.com/matlabcentral/answers/422244-why-do-buttons-apps-or-the-editor-in-matlab-respond-slowly-or-hang-on-macos>_`
 
