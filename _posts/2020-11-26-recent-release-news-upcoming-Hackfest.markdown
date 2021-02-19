---
layout: post
title:  "Recent release, news & upcoming Hackfest"
date:   2020-11-26 15:57:26 +0100
categories: jekyll update
author: Niklas
---

We are happy to announce our second major release: [**matRad "Blaise" (v2.10.1)**](https://github.com/e0404/matRad/releases/tag/v2.10.1).    
But before we swamp you with information about what’s changed and our future plans (which you can find further below), we want to direct your attention to something imminent:

**The upcoming (virtual) matRad Hackfest:**     
This yearly matRad coding session will take place next week on December 3 as virtual online hackathon. From 10:00 am (CET) on our group at DKFZ will dedicate this day to matRad development and give you the opportunity to join us via Zoom. You are free to use this opportunity to discuss projects, contributions, ask for help, request features, etc. You will be able to access the zoom room using (URL removed) on that day. You can also let us know beforehand to make sure we are not on a break when you join. We will be available to approx. 10:00 pm (CET).

And now for the remaining release information and news…

**Changes in the development team:**    
Previous matRad chief developer Mark passed on the baton to Niklas. While Mark is still involved in the development of matRad and you will occasionally find him contributing and replying to inquiries on GitHub, Niklas will now be the go-to person if you need help with your research projects in matRad.

**Release Info for matRad “Blaise”:**   
Download:    
[https://github.com/e0404/matRad/releases/tag/v2.10.1](https://github.com/e0404/matRad/releases/tag/v2.10.1)     
matRad 2.10.1 already marks the first patch to the original “Blaise” release 2.10.0 (which we allowed to breath on GitHub for some time in solitude), introducing quite a lot of changes, bug fixes and, primarily, the following new main features:    
- _Two open source Monte Carlo dose calculation engines have been integrated:_   
You can now call ompMC for photons (developed at the Pontifical Catholic University of Chile), and MCsquare for protons ([http://openmcsquare.org/](http://openmcsquare.org/), developed at Université catholique de Louvain)    
- _Independent dose and CT resolutions:_       
You can now define a resolution specific to the dose calculation and optimization accuracy to explore accuracy-to-runtime trade-offs.     
- _New optimization interface:_        
The optimization has been refactored to an object-oriented design with allows easy integration of your own objective & constraint functions, optimization problem structures and optimizers. Besides IPOPT, optimization is now also possible with the "fmincon" interior-point method from Matlab's optimization toolbox.    
- _DICOM exporter:_     
You can now export CT slices, RTstruct and RTDose from matRad (requires the image processing toolbox) using the new [```matRad_DicomExporter```](https://github.com/e0404/matRad/tree/v2.10.1/dicom/%40matRad_DicomExporter)    
- _Continuous Integration with TravisCI and Azure DevOps:_       
Changes and external contributions (i.e. Pull Requests) to matRad on github are now automatically tested in Matlab and Octave. This will ensure more stability of the code in the future.     

For a full list of changes see the [```ChangeLog.txt```](https://github.com/e0404/matRad/blob/v2.10.1/ChangeLog.txt) file in the matRad Code.

**What about future plans?**      
Currently, we are working towards version 3 which will include some of the following:       
- A multi-scenario framework for uncertainty quantification, robust optimization and 4D dose calculation & optimization.    
- Helium (physical & biological) dose calculation    
- Advanced biological modeling (including variable RBE proton models and a Helium model)     
- A new, modular and extended graphical user interface    
- Many more functionalities for Monte Carlo dose calculation (including an interface for TOPAS users!) 
If you are already very interested and want to try these functionalities out, let us know and we can point you to the right development branches.


P.S.:    
Not yet on our user-map [http://map.matrad.org](http://map.matrad.org)? Just let us know and we will add you asap!  
