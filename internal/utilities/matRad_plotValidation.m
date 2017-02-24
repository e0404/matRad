load('\\Mac\Home\Documents\Heidelberg\matRad validation\photons\matRadValPhanPhotons\0.mat');


Test = matRad_calcDoseDirect(ct,stf,pln,cst,1.3067);

resultGUI.refDose = Test.physicalDose;

matRadGUI