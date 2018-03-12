
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
```matlab
load 'LIVER.mat'

pln.bixelWidth      = 5;     % [mm] 
pln.gantryAngles    = [300]; % [�]
pln.couchAngles     = [0];   % [�]
pln.radiationMode   = 'carbon';
pln.bioOptimization = 'LEMIV_RBExD';

stf = matRad_generateStf(ct,cst,pln);

dij = matRad_calcParticleDose(ct,stf,pln,cst);

resultGUI = matRad_fluenceOptimization(dij,cst,pln);

matRadGUI
```
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

## More Code Examples on
## <span style="color:rgb(0,107,182)">https://github.com/e0404/matRad/tree/master/examples</span> 
---