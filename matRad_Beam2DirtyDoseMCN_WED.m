%% Testing the function
clear; close all; clc;

%% Proton Optimization 2 beam without dirty Dose
matRad_rc;
clear("dij","pln","resultGUI","ct","cst","stf")
load("PROSTATE.mat")

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'MCN'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIRef = resultGUI;

save("Reference_2BeamMCN.mat","resultGUI","pln","dij","ct","cst",'-mat')
totalphysDose = sum(resultGUIRef.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIRef.physicalDose);
Ref = matRad_calcQualityIndicators(cst,pln,resultGUIRef.dirtyDose);
sumDirtyDose = sum(resultGUIRef.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIRef.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Reference_2BeamMCN.mat","totalphysDose","physDose","Ref","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Bladder overdosing dirty Dose with penalty 100
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'MCN'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver100 = resultGUI;

save("BladderOver100_2BeamMCN.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIOver100.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver100.physicalDose);
Over100 = matRad_calcQualityIndicators(cst,pln,resultGUIOver100.dirtyDose);
sumDirtyDose = sum(resultGUIOver100.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver100.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BladderOver100_2BeamMCN.mat","totalphysDose","physDose","Over100","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Bladder overdosing dirty Dose with penalty 300
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(300,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'MCN'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver300 = resultGUI;

save("BladderOver300_2BeamMCN.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIOver300.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver300.physicalDose);
Over300 = matRad_calcQualityIndicators(cst,pln,resultGUIOver300.dirtyDose);
sumDirtyDose = sum(resultGUIOver300.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver300.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BladderOver300_2BeamMCN.mat","totalphysDose","physDose","Over300","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Bladder overdosing dirty Dose with penalty 500
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(500,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'MCN'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver500 = resultGUI;

save("BladderOver500_2BeamMCN.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIOver500.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver500.physicalDose);
Over500 = matRad_calcQualityIndicators(cst,pln,resultGUIOver500.dirtyDose);
sumDirtyDose = sum(resultGUIOver500.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver500.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BladderOver500_2BeamMCN.mat","totalphysDose","physDose","Over500","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Bladder overdosing dirty Dose with penalty 700
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(700,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'MCN'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver700 = resultGUI;

save("BladderOver700_2BeamMCN.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIOver700.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver700.physicalDose);
Over700 = matRad_calcQualityIndicators(cst,pln,resultGUIOver700.dirtyDose);
sumDirtyDose = sum(resultGUIOver700.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver700.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BladderOver700_2BeamMCN.mat","totalphysDose","physDose","Over700","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Bladder overdosing dirty Dose with penalty 900
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(900,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'MCN'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver900 = resultGUI;

save("BladderOver900_2BeamMCN.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIOver900.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver900.physicalDose);
Over900 = matRad_calcQualityIndicators(cst,pln,resultGUIOver900.dirtyDose);
sumDirtyDose = sum(resultGUIOver900.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver900.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BladderOver900_2BeamMCN.mat","totalphysDose","physDose","Over900","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Target underdosing dirty Dose with penalty 100 dmax 24
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(100,24));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'MCN'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIUnder100 = resultGUI;

save("TargetUnder100_2BeamMCN.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIUnder100.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIUnder100.physicalDose);
Under100 = matRad_calcQualityIndicators(cst,pln,resultGUIUnder100.dirtyDose);
sumDirtyDose = sum(resultGUIUnder100.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIUnder100.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("TargetUnder100_2BeamMCN.mat","totalphysDose","physDose","Under100","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Target underdosing dirty Dose with penalty 300 dmax 24
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(300,24));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'MCN'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIUnder300 = resultGUI;

save("TargetUnder300_2BeamMCN.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIUnder300.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIUnder300.physicalDose);
Under300 = matRad_calcQualityIndicators(cst,pln,resultGUIUnder300.dirtyDose);
sumDirtyDose = sum(resultGUIUnder300.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIUnder300.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("TargetUnder300_2BeamMCN.mat","totalphysDose","physDose","Under300","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Target underdosing dirty Dose with penalty 500 dmax 24
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(500,24));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'MCN'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIUnder500 = resultGUI;

save("TargetUnder500_2BeamMCN.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIUnder500.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIUnder500.physicalDose);
Under500 = matRad_calcQualityIndicators(cst,pln,resultGUIUnder500.dirtyDose);
sumDirtyDose = sum(resultGUIUnder500.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIUnder500.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("TargetUnder500_2BeamMCN.mat","totalphysDose","physDose","Under500","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Target underdosing dirty Dose with penalty 700 dmax 24
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(700,24));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'MCN'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIUnder700 = resultGUI;

save("TargetUnder700_2BeamMCN.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIUnder700.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIUnder700.physicalDose);
Under700 = matRad_calcQualityIndicators(cst,pln,resultGUIUnder700.dirtyDose);
sumDirtyDose = sum(resultGUIUnder700.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIUnder700.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("TargetUnder700_2BeamMCN.mat","totalphysDose","physDose","Under700","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Target underdosing dirty Dose with penalty 900 dmax 24
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(900,24));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'MCN'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIUnder900 = resultGUI;

save("TargetUnder900_2BeamMCN.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIUnder900.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIUnder900.physicalDose);
Under900 = matRad_calcQualityIndicators(cst,pln,resultGUIUnder900.dirtyDose);
sumDirtyDose = sum(resultGUIUnder900.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIUnder900.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("TargetUnder900_2BeamMCN.mat","totalphysDose","physDose","Under900","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Body mean dirty Dose with penalty 100 dmax 20
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{9,6}{2} = struct(DirtyDoseObjectives.matRad_MeanDirtyDose(100,20));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'MCN'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIMean100 = resultGUI;

save("BodyMean100_2BeamMCN.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIMean100.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIMean100.physicalDose);
Mean100 = matRad_calcQualityIndicators(cst,pln,resultGUIMean100.dirtyDose);
sumDirtyDose = sum(resultGUIMean100.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIMean100.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BodyMean100_2BeamMCN.mat","totalphysDose","physDose","Mean100","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Body mean dirty Dose with penalty 300 dmax 20
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{9,6}{2} = struct(DirtyDoseObjectives.matRad_MeanDirtyDose(300,20));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'MCN'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIMean300 = resultGUI;

save("BodyMean300_2BeamMCN.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIMean300.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIMean300.physicalDose);
Mean300 = matRad_calcQualityIndicators(cst,pln,resultGUIMean300.dirtyDose);
sumDirtyDose = sum(resultGUIMean300.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIMean300.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BodyMean300_2BeamMCN.mat","totalphysDose","physDose","Mean300","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Body mean dirty Dose with penalty 500 dmax 20
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{9,6}{2} = struct(DirtyDoseObjectives.matRad_MeanDirtyDose(500,20));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'MCN'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIMean500 = resultGUI;

save("BodyMean500_2BeamMCN.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIMean500.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIMean500.physicalDose);
Mean500 = matRad_calcQualityIndicators(cst,pln,resultGUIMean500.dirtyDose);
sumDirtyDose = sum(resultGUIMean500.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIMean500.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BodyMean500_2BeamMCN.mat","totalphysDose","physDose","Mean500","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Body mean dirty Dose with penalty 700 dmax 20
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{9,6}{2} = struct(DirtyDoseObjectives.matRad_MeanDirtyDose(700,20));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'MCN'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIMean700 = resultGUI;

save("BodyMean700_2BeamMCN.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIMean700.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIMean700.physicalDose);
Mean700 = matRad_calcQualityIndicators(cst,pln,resultGUIMean700.dirtyDose);
sumDirtyDose = sum(resultGUIMean700.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIMean700.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BodyMean700_2BeamMCN.mat","totalphysDose","physDose","Mean700","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Body mean dirty Dose with penalty 900 dmax 20
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{9,6}{2} = struct(DirtyDoseObjectives.matRad_MeanDirtyDose(900,20));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'MCN'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIMean900 = resultGUI;

save("BodyMean900_2BeamMCN.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIMean900.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIMean900.physicalDose);
Mean900 = matRad_calcQualityIndicators(cst,pln,resultGUIMean900.dirtyDose);
sumDirtyDose = sum(resultGUIMean900.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIMean900.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BodyMean900_2BeamMCN.mat","totalphysDose","physDose","Mean900","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Testing the function
clear; close all; clc;

%% Proton Optimization 2 beam without dirty Dose
matRad_rc;
clear("dij","pln","resultGUI","ct","cst","stf")
load("PROSTATE.mat")

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'WED'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIRef = resultGUI;

save("Reference_2BeamWED.mat","resultGUI","pln","dij","ct","cst",'-mat')
totalphysDose = sum(resultGUIRef.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIRef.physicalDose);
Ref = matRad_calcQualityIndicators(cst,pln,resultGUIRef.dirtyDose);
sumDirtyDose = sum(resultGUIRef.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIRef.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("Reference_2BeamWED.mat","totalphysDose","physDose","Ref","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Bladder overdosing dirty Dose with penalty 100
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'WED'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver100 = resultGUI;

save("BladderOver100_2BeamWED.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIOver100.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver100.physicalDose);
Over100 = matRad_calcQualityIndicators(cst,pln,resultGUIOver100.dirtyDose);
sumDirtyDose = sum(resultGUIOver100.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver100.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BladderOver100_2BeamWED.mat","totalphysDose","physDose","Over100","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Bladder overdosing dirty Dose with penalty 300
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(300,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'WED'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver300 = resultGUI;

save("BladderOver300_2BeamWED.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIOver300.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver300.physicalDose);
Over300 = matRad_calcQualityIndicators(cst,pln,resultGUIOver300.dirtyDose);
sumDirtyDose = sum(resultGUIOver300.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver300.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BladderOver300_2BeamWED.mat","totalphysDose","physDose","Over300","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Bladder overdosing dirty Dose with penalty 500
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(500,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'WED'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver500 = resultGUI;

save("BladderOver500_2BeamWED.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIOver500.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver500.physicalDose);
Over500 = matRad_calcQualityIndicators(cst,pln,resultGUIOver500.dirtyDose);
sumDirtyDose = sum(resultGUIOver500.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver500.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BladderOver500_2BeamWED.mat","totalphysDose","physDose","Over500","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Bladder overdosing dirty Dose with penalty 700
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(700,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'WED'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver700 = resultGUI;

save("BladderOver700_2BeamWED.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIOver700.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver700.physicalDose);
Over700 = matRad_calcQualityIndicators(cst,pln,resultGUIOver700.dirtyDose);
sumDirtyDose = sum(resultGUIOver700.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver700.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BladderOver700_2BeamWED.mat","totalphysDose","physDose","Over700","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Bladder overdosing dirty Dose with penalty 900
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(900,0));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'WED'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIOver900 = resultGUI;

save("BladderOver900_2BeamWED.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIOver900.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIOver900.physicalDose);
Over900 = matRad_calcQualityIndicators(cst,pln,resultGUIOver900.dirtyDose);
sumDirtyDose = sum(resultGUIOver900.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIOver900.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BladderOver900_2BeamWED.mat","totalphysDose","physDose","Over900","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Target underdosing dirty Dose with penalty 100 dmax 24
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(100,24));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'WED'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIUnder100 = resultGUI;

save("TargetUnder100_2BeamWED.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIUnder100.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIUnder100.physicalDose);
Under100 = matRad_calcQualityIndicators(cst,pln,resultGUIUnder100.dirtyDose);
sumDirtyDose = sum(resultGUIUnder100.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIUnder100.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("TargetUnder100_2BeamWED.mat","totalphysDose","physDose","Under100","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Target underdosing dirty Dose with penalty 300 dmax 24
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(300,24));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'WED'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIUnder300 = resultGUI;

save("TargetUnder300_2BeamWED.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIUnder300.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIUnder300.physicalDose);
Under300 = matRad_calcQualityIndicators(cst,pln,resultGUIUnder300.dirtyDose);
sumDirtyDose = sum(resultGUIUnder300.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIUnder300.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("TargetUnder300_2BeamWED.mat","totalphysDose","physDose","Under300","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Target underdosing dirty Dose with penalty 500 dmax 24
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(500,24));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'WED'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIUnder500 = resultGUI;

save("TargetUnder500_2BeamWED.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIUnder500.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIUnder500.physicalDose);
Under500 = matRad_calcQualityIndicators(cst,pln,resultGUIUnder500.dirtyDose);
sumDirtyDose = sum(resultGUIUnder500.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIUnder500.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("TargetUnder500_2BeamWED.mat","totalphysDose","physDose","Under500","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Target underdosing dirty Dose with penalty 700 dmax 24
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(700,24));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'WED'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIUnder700 = resultGUI;

save("TargetUnder700_2BeamWED.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIUnder700.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIUnder700.physicalDose);
Under700 = matRad_calcQualityIndicators(cst,pln,resultGUIUnder700.dirtyDose);
sumDirtyDose = sum(resultGUIUnder700.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIUnder700.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("TargetUnder700_2BeamWED.mat","totalphysDose","physDose","Under700","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Target underdosing dirty Dose with penalty 900 dmax 24
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(900,24));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'WED'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIUnder900 = resultGUI;

save("TargetUnder900_2BeamWED.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIUnder900.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIUnder900.physicalDose);
Under900 = matRad_calcQualityIndicators(cst,pln,resultGUIUnder900.dirtyDose);
sumDirtyDose = sum(resultGUIUnder900.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIUnder900.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("TargetUnder900_2BeamWED.mat","totalphysDose","physDose","Under900","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Body mean dirty Dose with penalty 100 dmax 20
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{9,6}{2} = struct(DirtyDoseObjectives.matRad_MeanDirtyDose(100,20));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'WED'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIMean100 = resultGUI;

save("BodyMean100_2BeamWED.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIMean100.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIMean100.physicalDose);
Mean100 = matRad_calcQualityIndicators(cst,pln,resultGUIMean100.dirtyDose);
sumDirtyDose = sum(resultGUIMean100.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIMean100.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BodyMean100_2BeamWED.mat","totalphysDose","physDose","Mean100","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Body mean dirty Dose with penalty 300 dmax 20
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{9,6}{2} = struct(DirtyDoseObjectives.matRad_MeanDirtyDose(300,20));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'WED'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIMean300 = resultGUI;

save("BodyMean300_2BeamWED.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIMean300.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIMean300.physicalDose);
Mean300 = matRad_calcQualityIndicators(cst,pln,resultGUIMean300.dirtyDose);
sumDirtyDose = sum(resultGUIMean300.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIMean300.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BodyMean300_2BeamWED.mat","totalphysDose","physDose","Mean300","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Body mean dirty Dose with penalty 500 dmax 20
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{9,6}{2} = struct(DirtyDoseObjectives.matRad_MeanDirtyDose(500,20));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'WED'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIMean500 = resultGUI;

save("BodyMean500_2BeamWED.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIMean500.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIMean500.physicalDose);
Mean500 = matRad_calcQualityIndicators(cst,pln,resultGUIMean500.dirtyDose);
sumDirtyDose = sum(resultGUIMean500.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIMean500.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BodyMean500_2BeamWED.mat","totalphysDose","physDose","Mean500","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Body mean dirty Dose with penalty 700 dmax 20
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{9,6}{2} = struct(DirtyDoseObjectives.matRad_MeanDirtyDose(700,20));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'WED'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIMean700 = resultGUI;

save("BodyMean700_2BeamWED.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIMean700.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIMean700.physicalDose);
Mean700 = matRad_calcQualityIndicators(cst,pln,resultGUIMean700.dirtyDose);
sumDirtyDose = sum(resultGUIMean700.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIMean700.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BodyMean700_2BeamWED.mat","totalphysDose","physDose","Mean700","sumDirtyDose","maxDirtyDose",'-mat','-append')

%% Proton Optimization 2 beam Body mean dirty Dose with penalty 900 dmax 20
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{9,6}{2} = struct(DirtyDoseObjectives.matRad_MeanDirtyDose(900,20));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% Define the flavor of optimization
quantityOpt   = 'RBExD'; 
modelName     = 'WED'; %MCN for protons, HEL for helium, LEM for carbon

pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,'nomScen');

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUIMean900 = resultGUI;

save("BodyMean900_2BeamWED.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUIMean900.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUIMean900.physicalDose);
Mean900 = matRad_calcQualityIndicators(cst,pln,resultGUIMean900.dirtyDose);
sumDirtyDose = sum(resultGUIMean900.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUIMean900.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BodyMean900_2BeamWED.mat","totalphysDose","physDose","Mean900","sumDirtyDose","maxDirtyDose",'-mat','-append')

