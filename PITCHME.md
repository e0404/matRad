## matRad 

![Image](https://github.com/e0404/matRad/blob/master/dicomImport/matrad_logo.png)
An open source software for radiation treatment planning of intensity-modulated photon, proton, and carbon ion therapy.

---

### matRad provides functionalites for 

- Ray tracing
- Photon dose calculation
- Proton dose calculation
- Carbon ion dose calculation 
- Inverse planning 
- Multileaf collimator sequencing
- Basic treatment plan visualization and evaluation

---
### Code

```
stf       = matRad_generateStf(ct,cst,pln);
dij       = matRad_calcParticleDose(ct,stf,pln,cst);
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
matRadGUI
```
