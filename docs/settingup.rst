.. toctree::
   :maxdepth: 2
   :hidden:

   guide/matradscript
   guide/gui

.. |matRad_logo_header| image:: ../matRad/gfx/matRad_logo.png
   :width: 130 px
   :alt: matRad
   :target: https://www.matRad.org

.. |matRad_logo| image:: ../matRad/gfx/matRad_logo.png
   :width: 60 px
   :alt: matRad
   :target: https://www.matRad.org

.. _settingup:

===============================
Setting Up |matRad_logo_header|
===============================

To get |matRad_logo| running you have two options:

1. Source Code for people with Matlab installation:
---------------------------------------------------

If you have MATLAB you can just `get a local copy of the source code <How-to-get-a-local-copy>`_. Then you have to choose whether you want to `use the  GUI <How-to-execute-matRadGUI>`_ or `execute the main script <How-to-execute-matRad>`_.

2. Standalone installation:
---------------------------
If you do not have MATLAB installed on your PC you are restricted to use |matRad_logo|'s standalone application. For this option a MATLAB installation is not required. The installer is available for Windows, Mac and Linux `with the latest release <https://github.com/e0404/matRad/releases/latest>`_.

Steps for installation:
1. Download the |matRad_logo| installer of the latest `Release <https://github.com/e0404/matRad/releases/latest>`_ to your system (`Win64 <https://github.com/e0404/matRad/releases/download/v2.10.1/matRad_installerWin64_v2.10.1.exe>`_, `Mac64 <https://github.com/e0404/matRad/releases/download/v2.10.1/matRad_installerMac64_v2.10.1.dmg>`_, `Linux64 <https://github.com/e0404/matRad/releases/download/v2.10.1/matRad_installerLinux64_v2.10.1.install>`_)

2. Run the respective installer for your system

    * **Windows**: Run the downloaded executable installer
    * **Linux**: Run the executable install script. Make sure that the ``*.install`` file has executable permissions.
    * **Mac**: Here we provide a dmg containing the installer (Since the installer is not Apple-certified, you might explicitly launch it from the terminal or by right-click and then open).

    After that, you should be guided through the installation process.
    Note that the installers will want to download the "Matlab Runtime" from Mathworks in the process. The runtime is quite large (~2GB) and is required to run compiled deployed applications written in Matlab.

3. Run |matRad_logo|:

    * **Windows**: Just like with every other program, you should have a desktop icon.
    * **Mac**: Per default |matRad_logo| will be installed to ``/Applications/matRad``. To run |matRad_logo| navigate to ``/Applications/matRad/applications`` and double click or right click -> open on the ``matRad.app`` application. The first startup might take a few seconds.

      .. note::
          An installation warning appears that |matRad_logo| is from an unverified developer. You can solve this issue by opening the installer from the context menu (depending on the configuration either Ctrl + click or right-click on the icon, and then click "Open" in the menu). You will then get the option to open the application in spite of the missing verification and thus to install |matRad_logo|.

    * **Linux & Mac**: To start |matRad_logo|, you can alternatively use the provided ``run_matRad.sh`` script from the terminal. It requires one argument which gives the path to the installed Matlab-Runtime. Refer to the ``readme_linux.txt`` and ``readme_mac.txt`` in your installation directory for more information.

Patient/Phantom files
---------------------
The patient files should be included with the installer and will be installed into the desired location. For Windows, for example, they can be found within the "application" folder of the chosen installation directory. Moreover, we also provide extra links to the open source patient files stored in |matRad_logo|'s native ``*.mat`` format for a `head and neck case <https://github.com/e0404/matRad/blob/master/phantoms/HEAD_AND_NECK.mat>`_, a `liver case <https://github.com/e0404/matRad/blob/master/phantoms/LIVER.mat>`_, a `prostate case <https://github.com/e0404/matRad/blob/master/phantoms/PROSTATE.mat>`_, a `box phantom <https://github.com/e0404/matRad/blob/master/phantoms/BOXPHANTOM.mat>`_, and AAPM's `TG119 phantom <https://github.com/e0404/matRad/blob/master/phantoms/TG119.mat>`_. Otherwise you need to start with a `DICOM import of your own patient data <The-dicom-import>`_.


An overview of the functions of the |matRad_logo| GUI can be found `here <The-matRad-GUI>`_.

If you want to import patient data, check out `the dicom import functionality <The-dicom-import>`_.