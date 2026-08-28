.. _result:

============================
The resultGUI Data Structure
============================

The optimization output is declared as ``resultGUI`` and contains at least two fields. The first field ``w`` contains the optimized weights, and the second field, ``physicalDose``, holds the physical dose cube in CT dimensions :math:`N_x \times N_y \times N_z` (:math:`N_x, N_y, N_z` = number of voxels in x-, y-, and z-direction).
Each field in the ``resultGUI`` struct can be easily accessed via MATLAB's dot-operator. For instance, ``resultGUI.w`` outputs the fluence weight vector.

In case of performing a biological optimization using carbon ions, the ``resultGUI`` struct holds several additional cubes as fields:

- ``resultGUI.effect`` - represents the biological effect cube
- ``resultGUI.RBExDose`` - contains the biological effective dose = RBE × physical dose
- ``resultGUI.RBE`` - holds the relative biological effectiveness cube
- ``resultGUI.alpha`` - represents the radiosensitivity parameter :math:`\alpha_\mathrm{Particle}` cube from the linear quadratic model
- ``resultGUI.beta``  - represents the radiosensitivity parameter :math:`\beta_\mathrm{Particle}` cube from the linear quadratic model

Executing ``matRadGUI`` in MATLAB's command window will start the GUI and display all available data of ``resultGUI``.
