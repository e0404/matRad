.. _cort:

================
The CORT Dataset
================

The matRad release includes five image datasets of phantoms and patients. 
Consequently, everybody can start treatment planning out of the box without the need for their own data. 
Besides a simple cubic box phantom, matRad also features the four datasets that have been published as the `CORT dataset <http://dx.doi.org/10.1186/2047-217X-3-37>`_.

The `CORT dataset <http://dx.doi.org/10.1186/2047-217X-3-37>`_ (common optimization for radiation therapy) is an open dataset intended to be used by researchers when developing and contrasting radiation treatment planning optimization algorithms. It is comprised of datasets for a prostate case, a liver case, a head and neck case, and a standard IMRT phantom (TG119). In matRad, these datasets are already imported and ready to be loaded as MATLAB native ``.mat`` files.
To allow rapid prototyping and reasonable performance in educational settings, the CT images were downsampled at import. 

If you wish to take a look at the original data, you can find it `here <http://dx.doi.org/10.5524/100110>`_.  
The dataset contains the original Digital Imaging and Communications in Medicine (DICOM) computed tomography (CT) scan, as well as the DICOM structure file. In addition, the dose-influence matrix from a variety of beam/couch angle pairs is provided for each case.

Besides the pre-configured patient datasets, matRad also features a :ref:`DICOM import <dicom>`.