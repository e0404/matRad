%% Dose depending on LET
clear all; close all; clc;

%% Proton Optimization with one beam

matRad_rc;
load('PROSTATE.mat')

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
% pln.propStf.gantryAngles  = [90 270];
pln.propStf.gantryAngles  = 90;
% pln.propStf.couchAngles   = [0 0];
pln.propStf.couchAngles   = 0;
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

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);

%% User Interface
matRadGUI;

%% Searching for Voxel
% Info: 
% Cube Index: [89 66 42 ]
% Space Coordinates: [198 267 126]
% HU Value: 58.6

%numVoxel = ct.cube{1}(89,66,42);
% roundnumVoxel = round(numVoxel);
ind = sub2ind(size(ct.cube{1}),89,66,42);

%d = 1;
%while dij.mLETDose{1}(d,1) == 0
%    d = d+1;
%end

rowLET = full(dij.mLETDose{1}(ind,:));
rowPhysDose = full(dij.physicalDose{1}(ind,:));

% Coordinates = ct.cube{1}(198,267,126)
% value = ct.x(1,ct.x==198)
% row = dij.mLETDose(dij.mLETDose{1} == voxel,1)
% [value,position] = nonzeros(dij.mLETDose{1})

%% mLETdose

mLETdose = rowLET.*resultGUI.w';
rowPhysicalDose = rowPhysDose.*resultGUI.w';

% mLETdose(mLETdose==0) = [];

% plot LETd depending on bixels j
% bixel = row./numVoxel;
% j = bixel';
figure
plot(mLETdose)
numnon0LET = nnz(mLETdose); 
mLETdose(mLETdose==0)=[];
sizeLETdose = size(mLETdose);
rowPhysicalDose(rowPhysicalDose==0) = [];

LET = mLETdose./rowPhysicalDose;

LET_thres = 2.5;

dirtydose = LET(LET>LET_thres);
cleandose = LET(LET<LET_thres);

figure
h = histogram(LET);
xlabel('LET'); ylabel('bixel')
Nbins = morebins(h);
h.NumBins = length(LET);

% sortphysicalDose = sort(rowPhysDose);

% figure
% h1 = histogram(rowPhysDose,LET);
% 
% figure
% h2 = histogram2(LET,rowPhysDose);
% Nbins = morebins(h2);
% h2.NumBins = length(LET);

% figure
% scatter(mLETdose,rowPhysDose)
%histogram(mLETdose,rowPhysDose,)
% only1weight = mLETdose(1,:);
% figure
% plot(only1weight)

%% LET depending on Dose

% ixProfileY = round(pln.propStf.isoCenter(1,2)./ct.resolution.y);
% slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
% doseslice = resultGUI.physicalDose(:,ixProfileY,slice);
% 
% figure
% plot(mLETdose(1:183,1),doseslice)

% histogram
figure
histogram(mLETdose)

%% A different voxel row

% while dij.mLETDose{1}(d,1) ~= 0
%     d = d+1;
% end

% Info:
% Cube Index: [92 74 42]

ind2 = sub2ind(size(ct.cube{1}),92,74,42);
rowLET2 = full(dij.mLETDose{1}(ind2,:));

mLETdose2 = rowLET2.*resultGUI.w';

figure
plot(mLETdose2)

% histogram
figure
histogram(mLETdose2)

%% A different voxel row

% Info:
% Cube Index: [92 98 42]

ind3 = sub2ind(size(ct.cube{1}),92,98,42);
rowLET3 = full(dij.mLETDose{1}(ind3,:));

mLETdose3 = rowLET3.*resultGUI.w';

figure
plot(mLETdose3)

% histogram
figure
histogram(mLETdose3)


%% Proton Optimization with two beams

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
%pln.propStf.gantryAngles  = 90;
pln.propStf.couchAngles   = [0 0];
%pln.propStf.couchAngles   = 0;
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
dij= matRad_calcParticleDose(ct,stf,pln,cst);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);

%% Choosing a row
% t = 1;
% while dij.mLETDose{1}(t,1) == 0
%     t = t+1;
% end

% Info:
% Cube Index: [108 71 42]

ind_d = sub2ind(size(ct.cube{1}),108,71,42);
rowLET_d = full(dij.mLETDose{1}(ind_d,:));

mLETdose_d = rowLET_d.*resultGUI.w';

figure
plot(mLETdose_d)

% histogram
figure
histogram(mLETdose_d)


%% A different voxel row
ind_d2 = sub2ind(size(ct.cube{1}),92,74,42);
rowLET_d2 = full(dij.mLETDose{1}(ind_d2,:));

mLETdose_d2 = rowLET_d2.*resultGUI.w';

figure
plot(mLETdose_d2)

% histogram
figure
histogram(mLETdose_d2)

%% A different voxel row
ind_d3 = sub2ind(size(ct.cube{1}),92,98,42);
rowLET_d3 = full(dij.mLETDose{1}(ind_d3,:));

mLETdose_d3 = rowLET_d3.*resultGUI.w';

figure
plot(mLETdose_d3)

% histogram
figure
histogram(mLETdose_d3)

%% Proton Optimization with three beams

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 120 270];
%pln.propStf.gantryAngles  = 90;
pln.propStf.couchAngles   = [0 0 0];
%pln.propStf.couchAngles   = 0;
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
dij= matRad_calcParticleDose(ct,stf,pln,cst);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);

%% Choosing a row
% t = 1;
% while dij.mLETDose{1}(t,1) == 0
%     t = t+1;
% end

% Info:
% Cube Index: [108 71 42]

ind_t = sub2ind(size(ct.cube{1}),108,71,42);
rowLET_t = full(dij.mLETDose{1}(ind_t,:));

mLETdose_t = rowLET_t.*resultGUI.w';

figure
plot(mLETdose_t)

% histogram
figure
histogram(mLETdose_t)


%% A different voxel row
ind_t2 = sub2ind(size(ct.cube{1}),92,74,42);
rowLET_t2 = full(dij.mLETDose{1}(ind_t2,:));

mLETdose_t2 = rowLET_t2.*resultGUI.w';

figure
plot(mLETdose_t2)

% histogram
figure
histogram(mLETdose_t2)

%% A different voxel row
ind_t3 = sub2ind(size(ct.cube{1}),92,98,42);
rowLET_t3 = full(dij.mLETDose{1}(ind_t3,:));

mLETdose_t3 = rowLET_t3.*resultGUI.w';

figure
plot(mLETdose_t3)

% histogram
figure
histogram(mLETdose_t3)

%% Proton Optimization with four beams

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 120 270 300];
%pln.propStf.gantryAngles  = 90;
pln.propStf.couchAngles   = [0 0 0 0];
%pln.propStf.couchAngles   = 0;
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
dij= matRad_calcParticleDose(ct,stf,pln,cst);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);

%% Choosing a row
% t = 1;
% while dij.mLETDose{1}(t,1) == 0
%     t = t+1;
% end

% Info:
% Cube Index: [108 71 42]

ind_f = sub2ind(size(ct.cube{1}),108,71,42);
rowLET_f = full(dij.mLETDose{1}(ind_f,:));

mLETdose_f = rowLET_f.*resultGUI.w';

figure
plot(mLETdose_f)

% histogram
figure
histogram(mLETdose_f)


%% A different voxel row
ind_f2 = sub2ind(size(ct.cube{1}),92,74,42);
rowLET_f2 = full(dij.mLETDose{1}(ind_f2,:));

mLETdose_f2 = rowLET_f2.*resultGUI.w';

figure
plot(mLETdose_f2)

% histogram
figure
histogram(mLETdose_f2)

%% A different voxel row
ind_f3 = sub2ind(size(ct.cube{1}),92,98,42);
rowLET_f3 = full(dij.mLETDose{1}(ind_f3,:));

mLETdose_f3 = rowLET_f3.*resultGUI.w';

figure
plot(mLETdose_f3)

% histogram
figure
histogram(mLETdose_f3)

%% Proton Optimization with five beams

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [60 90 120 270 300];
%pln.propStf.gantryAngles  = 90;
pln.propStf.couchAngles   = [0 0 0 0 0];
%pln.propStf.couchAngles   = 0;
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
dij= matRad_calcParticleDose(ct,stf,pln,cst);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);

%% Choosing a row
% t = 1;
% while dij.mLETDose{1}(t,1) == 0
%     t = t+1;
% end

% Info:
% Cube Index: [108 71 42]

ind_v = sub2ind(size(ct.cube{1}),108,71,42);
rowLET_v = full(dij.mLETDose{1}(ind_v,:));

mLETdose_v = rowLET_v.*resultGUI.w';

figure
plot(mLETdose_v)

% histogram
figure
histogram(mLETdose_v)


%% A different voxel row
ind_v2 = sub2ind(size(ct.cube{1}),92,74,42);
rowLET_v2 = full(dij.mLETDose{1}(ind_v2,:));

mLETdose_v2 = rowLET_v2.*resultGUI.w';

figure
plot(mLETdose_v2)

% histogram
figure
histogram(mLETdose_v2)

%% A different voxel row
ind_v3 = sub2ind(size(ct.cube{1}),92,98,42);
rowLET_v3 = full(dij.mLETDose{1}(ind_v3,:));

mLETdose_v3 = rowLET_v3.*resultGUI.w';

figure
plot(mLETdose_v3)

% histogram
figure
histogram(mLETdose_v3)

%% Proton Optimization with six beams

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [60 90 120 240 270 300];
%pln.propStf.gantryAngles  = 90;
pln.propStf.couchAngles   = [0 0 0 0 0 0];
%pln.propStf.couchAngles   = 0;
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
dij= matRad_calcParticleDose(ct,stf,pln,cst);

% Inverse Optimization for IMPT
resultGUI= matRad_fluenceOptimization(dij,cst,pln);

%% Choosing a row
% t = 1;
% while dij.mLETDose{1}(t,1) == 0
%     t = t+1;
% end

% Info:
% Cube Index: [108 71 42]

ind_s = sub2ind(size(ct.cube{1}),108,71,42);
rowLET_s = full(dij.mLETDose{1}(ind_s,:));

mLETdose_s = rowLET_s.*resultGUI.w';

figure
plot(mLETdose_s)

% histogram
figure
histogram(mLETdose_s)


%% A different voxel row
ind_s2 = sub2ind(size(ct.cube{1}),92,74,42);
rowLET_s2 = full(dij.mLETDose{1}(ind_s2,:));

mLETdose_s2 = rowLET_s2.*resultGUI.w';

figure
plot(mLETdose_s2)

% histogram
figure
histogram(mLETdose_s2)

%% A different voxel row
ind_s3 = sub2ind(size(ct.cube{1}),92,98,42);
rowLET_s3 = full(dij.mLETDose{1}(ind_s3,:));

mLETdose_s3 = rowLET_s3.*resultGUI.w';

figure
plot(mLETdose_s3)

% histogram
figure
histogram(mLETdose_s3)