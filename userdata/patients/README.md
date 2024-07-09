# Folder for user-imported patients

This folder can be used to store user-imported patients (default import location for DICOM and binary imports). *.mat files containing matRad structures (ct,cst,stf, etc.) stored in here will be added on to the search path once `matRad_rc` is run, but will by default be ignored within git.

For maximum compatibility, we suggest to store manually imported mat-files in the v7 format when using Matlab's `save` function.