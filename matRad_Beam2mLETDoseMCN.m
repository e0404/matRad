%% Proton Optimization 2 beam Bladder overdosing dirty Dose with penalty 100
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cst{8,6}{2} = struct(LETdObjectives.matRad_SquaredOverdosingLETd(100,0));

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
resultGUILETOver100 = resultGUI;

save("BladdermLETDoseOver100_2BeamMCN.mat","resultGUI","pln","dij","ct","cst",'-mat')

totalphysDose = sum(resultGUILETOver100.physicalDose(:));
physDose = matRad_calcQualityIndicators(cst,pln,resultGUILETOver100.physicalDose);
LETOver100 = matRad_calcQualityIndicators(cst,pln,resultGUILETOver100.dirtyDose);
sumDirtyDose = sum(resultGUILETOver100.dirtyDose(cst{8,4}{1}));
maxDirtyDose = max(resultGUILETOver100.dirtyDose(cst{8,4}{1})); % noch nicht berechnet
save("BladdermLETDoseOver100_2BeamMCN.mat","totalphysDose","physDose","LETOver100","sumDirtyDose","maxDirtyDose",'-mat','-append')