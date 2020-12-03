Readme for installation of the standalone on Windows

1. Install matRad

Two executables are available to install matRad on your machine 

• With Runtime

    >> matRad_installerWin64_wRT.exe

The installer for MATLAB Runtime version 9.7 (R2019b) is already
packaged in this installer.


• Without Runtime

    >> atRad_installerWin64.install

This installer will connect to the internet to download the installer for 
the MATLAB Runtime. For this, make sure that the terminal can connect to the 
internet (proxy settings are correctly set), and that the SSL CA certificate 
is installed on your machine (the installer might need root privileges to 
access the certificates). 

To run the installer, just execute the executable file.


The matRad Installer will then be launched. There you can 

    (1) Choose Destination folder for matRad

    (2) Choose Destination folder for MATLAB Runtime (if not yet installed)


If the environment is properly configured and the installer without Runtime 
still doesn't run, it will not be possible to install matRad using this 
executable. If you face this issue, use the installer with Runtime instead.


2. Launch matRad

To launch matRad, just start the executable matRad.exe from the installation directory as with any other program.
The patient files are located within the chosen installation directory in the "application" folder.



