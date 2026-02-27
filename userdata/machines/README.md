# Folder for user-configured machines

This folder can be used to store user-imported machines as *.mat files. Machines in here will be added on to the search path once `matRad_rc` is run and subsequently recognized in matRad including the GUI, but will by default be ignored within git.
For the format of machine files check-out the Generic machines in the matRad/basedata subfolder or consult the Wiki.

For maximum compatibility, we suggest to store manually imported mat-files in the v7 format when using Matlab's `save` function.