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
modelName     = 'constRBE'; %MCN for protons, HEL for helium, LEM for carbon

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
resultGUIRef= matRad_fluenceOptimization(dij,cst,pln);

%% Proton Optimization 2 beam Bladder overdosing dirty Dose with penalty 100
clear("dij","stf")

% instead of dirtyDose you can try LETxDose Objectives, just use LETxDose instead of DirtyDose inside the brackets, but for both "DirtyDose"
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(100,0));
 
% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUIOver100= matRad_fluenceOptimization(dij,cst,pln);

%% Proton Optimization 2 beam Bladder overdosing dirty Dose with penalty 300
clear("dij","stf")
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(300,0));

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUIOver300= matRad_fluenceOptimization(dij,cst,pln);

%% Proton Optimization 2 beam Bladder overdosing dirty Dose with penalty 500
clear("dij","stf")
cst{8,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredOverdosingDirtyDose(500,0));
% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUIOver500= matRad_fluenceOptimization(dij,cst,pln);

%% Proton Optimization 2 beam Target underdosing dirty Dose with penalty 100 dmax 24
clear("dij","cst","stf")

load("PROSTATE.mat")

cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(100,24));

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUIUnder100= matRad_fluenceOptimization(dij,cst,pln);

%% Proton Optimization 2 beam Target underdosing dirty Dose with penalty 300 dmax 24
clear("dij","stf")
cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(300,24));

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUIUnder300= matRad_fluenceOptimization(dij,cst,pln);

%% Proton Optimization 2 beam Target underdosing dirty Dose with penalty 500 dmax 24
clear("dij","stf")
cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(500,24));

% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);

% it only works with a Boxphantom and a LET threshold between 2.7 and 7.6:
dij = matRad_calcDirtyDose(2,dij);

% Inverse Optimization for IMPT
resultGUIUnder500 = matRad_fluenceOptimization(dij,cst,pln);

% plot distributions by using matRad_plotSliceWrapper
