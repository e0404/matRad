Readme for installation of the standalone on linux
1. Install matRad

Two executables are available to install matRad on your machine 

• With Runtime

    >> matRad_installerLinux64_wRT.install

The installer for MATLAB Runtime version 9.7 (R2019b) is already
packaged in this installer.


• Without Runtime

    >> matRad_installerLinux64.install

This installer will connect to the internet to download the installer for 
the MATLAB Runtime. For this, make sure that the terminal can connect to the 
internet (proxy settings are correctly set), and that the SSL CA certificate 
is installed on your machine (the installer might need root privileges to 
access the certificates). 

To run the installer from the terminal, go to the directory where the
installer is located, and type in the command prompt

       ./matRad_installerLinux64_wRT.install
or
       ./matRad_installerLinux64.install
or
       sudo ./matRad_installerLinux64.install


The matRad Installer will then be launched. There you can 

    (1) Choose Destination folder for matRad

    (2) Choose Destination folder for MATLAB Runtime


If the environment is properly configured and the installer without Runtime 
still doesn't run, it will not be possible to install matRad using this 
executable. If you face this issue, use the installer with Runtime instead.


Files and Folders created in the matRad directory after installation
====================================================================

- appdata

- application
    • example files
        >> BOXPHANTOM.mat
        >> HEAD_AND_NECK.mat
        >> LIVER.mat
        >> PROSTATE.mat
        >> TG119.mat

    • the shell script to launch matRad
        >> run_matRad.sh

- uninstall


2. Launch matRad

To run the shell script, type in the command prompt
   
       ./run_matRad.sh <mcr_directory> <argument_list>
       
    <mcr_directory> is the directory where version 9.7 of the  
    MATLAB Runtime is installed. 

    <argument_list> is all the arguments you want to pass to 
    your application. 

For example, if you have version 9.7 of the MATLAB Runtime installed 
in /mathworks/home/application/v97, run the shell script as
    
       ./run_matRad.sh /mathworks/home/application/v97

