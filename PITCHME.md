
![Logo](dicomImport/matrad_logo.png)

###  A radiation treatment planning software for intensity-modulated <span style="color:rgb(0,107,182); font-size: 1em;">photon</span>, <span style="color:rgb(0,107,182)">proton</span> and <span style="color:rgb(0,107,182)">carbon ion</span> therapy.
---
# <span style="color:rgb(0,107,182)">Free</span>  Software for reasearch and education
---
## matRad provides functionalites for 
- DICOM import
- ray tracing
- photon dose calculation
- proton & carbon dose calculation
- inverse planning 
- multileaf collimator sequencing
- treatment plan visualization and evaluation
---
### Graphical User Interface
![Logo](https://github.com/e0404/matRad/wiki/images/GUI-Guide_optimizedGUIScreenshot.png)
---
## Code Example
+++?gist=becker89/d1681e8ff3ba1e22dd26b645ad6b0544
@[1](import a open source liver patient)
@[3-7](define your treatment plan)
@[9](generate beam and ray geometry)
@[11](dose calculation - obtain dose influence matrix)
@[13](inverse planning for IMPT)
@[15](start GUI for visualization of result)
---
![Logo](https://github.com/e0404/matRad/wiki/images/matRadvalidation.png)
---
## Performance 
![Logo](https://github.com/e0404/matRad/wiki/images/matRadPerformanceTable.png)
---
## matRad webinar 
![Video](https://www.youtube.com/embed/40_n7BIqLdw)
---
## Get in touch
### <span style="color:rgb(0,107,182)">https://matRad.org/</span> 
### matRad@dkfz.de

---
## Code Example

```matlab
beamIx = 2;
stf(beamIx).isoCenter(1) = stf(beamIx).isoCenter(1) + 3;

resultGUI_isoShift = matRad_calcDoseDirect(ct,stf,pln,cst,w);

slice = round(stf.isoCenter(3)./ct.resolution.z);

DoseDiff = resultGUI(:,:,slice) - resultGUI_isoShift(:,:,slice);

figure, imagesc(DoseDiff);colorbar
```
@[1-2](Lets simulate a lateral displacement in x of the second beam)
@[4](recalculate dose using previously optimized beamlet weights)
@[6](determine axial iso center slice)
@[8](calculate dose difference)
@[10](plot dose difference slice)
---

```matlab
%% generate steering file
stf = matRad_generateStf(ct,cst,pln);

%% dose calculation
dij = matRad_calcPhotonDose(ct,stf,pln,cst);

%% inverse planning for imrt
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% indicator calculation
cst = matRad_indicatorWrapper(cst,pln,resultGUI);

%% sequencing
resultGUI = matRad_xiaLeafSequencing(resultGUI,stf,dij,5);

%% DAO
resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln);
matRad_visApertureInfo(resultGUI.apertureInfo);

%% start gui for visualization of result
matRadGUI
matRad_showDVH(cst,pln)
```
---
