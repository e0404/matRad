.. _dicomimport:

================
The DICOM import
================

Disclaimer
==========

Even though DICOM data is supposed to be highly standardized, it still features substantial variation. Consequently, it is very difficult for us to eliminate all bugs in matRad's DICOM import. As of June 2015, we have a beta version included in matRad which has been thoroughly tested using data from our center(s). We would like to remind all users to take utmost care when running matRad's DICOM import and double-check the imported CT data in the :ref:`ct-struct` and the imported structure sets in the :ref:`cst-cell`. If you do find any bugs, help us improve matRad's DICOM import and `drop us a line <https://github.com/e0404/matRad/issues>`_.

Importing patient data
======================

1. `Preparing patient data for import <#step1>`_
2. `Starting the import GUI <#step2>`_
3. `Completing the import <#step3>`_
4. `Checking the import <#step4>`_

.. _step1:

Step 1: Preparing patient data for import
-----------------------------------------

As of now, all ``*.dcm`` files need to be located in the same folder (no subfolders) for the import to work. That includes all CT images and the structure file(s). Dose series files enable you to directly import the corresponding dose distribution to your treatment plan. If you wish to include the treatment plan files and dose series files in your import, both need to be located in the same folder as well. It is possible to have DICOM data for multiple patients in this folder, as you are able to specify the patient, CT series, structure set, and treatment plan later during the import.

.. _step2:

Step 2: Starting the import GUI
-------------------------------

To start the import, you first need to add the matRad folder `dicom <https://github.com/e0404/matRad/tree/master/dicom>`_ to your path. You can do so either by *right-click* → *Add to Path* → *Selected Folders and Subfolders* or by typing ``addpath(genpath('dicom'))`` in your Command Window. ``genpath`` adds the subfolder *HLUT library* to your path, too.

Now you have access to the matRad functions designated for the DICOM import. To start the GUI you can simply type ``matRad_importDicomGUI`` in your Command Window. If you already work with matRad and have matRadGUI open, you can also start the import GUI with the *Load DICOM* button.

You should see something like this:

.. image:: /images/dicomImport/dicomImportGUI_empty.png

.. _step3:

Step 3: Completing the import
-----------------------------

Now you can choose your patient directory. Simply click the *Browse* button in the right upper corner and navigate to the desired directory. After choosing a folder, it will directly be analyzed for different patients.

The result will look something like this:

.. image:: /images/dicomImport/dicomImportGUI_withPatients.png

Now you can choose your patient, CT series, structure set, RT plan, and dose series. If there are several dose cubes referring to the same plan, you can decide to import just one, several, or all dose cubes.
In the bottom left corner, you can see the resolution of the chosen image series. You can adjust these values to import an interpolated cube with your specified resolution. Downsampling the CT makes sense in order to restrict matRad's memory usage. If your image series contains non-equally spaced slices, you have to specify a resolution. Another option is to use the resolution of the dose grid by activating the check box *Use RT Dose grid* after choosing a dose series.
Now you can click on *Import* and will see a progress bar pop-up.

Your command window will show you the following output:

.. image:: /images/dicomImport/dicomImportGUI_Output1.png
.. image:: /images/dicomImport/dicomImportGUI_Output2.png

If you import steering information from a proton or carbon treatment plan, you have to specify the base data of your machine as this information is not saved in DICOM files. Therefore, you will see a pop-up window in which you can navigate to and choose the correct base data for the beam quality you use. If you import steering information from a photon treatment plan, the generic photon base data will be selected automatically.

After the import has finished, a 'Save as'-dialogue will open to save the imported data in a ``*.mat`` file.

.. image:: /images/dicomImport/dicomImportGUI_savePatient.png

.. _step4:

Step 4: Checking the import
---------------------------

Objectives and constraints are not imported from DICOM files but default values will be assigned to them. If any structure is marked as *'tv'*, *'target'*, *'gtv'*, *'ctv'*, *'ptv'*, *'boost'*, or *'tumor'*, its VOI type will be set to *'TARGET'*, its priority to 1, its objective to *square deviation* with a penalty of 800 and a dose of 30 Gy. All other structures will be set to organs at risk (*OAR*) with a priority of 2.
The default values for biological planning are set to ``alphaX = 0.1`` and ``betaX = 0.05`` for all structures.
After importing and saving the imported files, you should check and adapt all objectives and constraints for your application.

If you have imported a photon treatment plan, you can now adjust the base data to your own machine by adapting :ref:`pln.machine`.

Now you can use this patient data set just like the ones provided with matRad.