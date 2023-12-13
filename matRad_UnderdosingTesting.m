%% Proton Optimization 2 beam without dirty Dose
matRad_rc;
clear("dij","pln","resultGUI","ct","cst","stf")
load("PROSTATE.mat")

cube = zeros(183,183,90);
cube(cst{6,4}{1}) = 1;
vResolution = ct.resolution;
vMargin = [];
vMargin.x = 5;
vMargin.y = 5;
vMargin.z = 5;
mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin);

cst{11,1}    = 10;
cst{11,2}    = 'Margin';
cst{11,3}    = 'OAR';
cst{11,4}{1} = find(mVOIEnlarged);
cst{11,5}    = cst{3,5};

cst{11,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,60)); 

cst{6,6}{2} = struct(mLETDoseObjectives.matRad_SquaredUnderdosingmLETDose(100,20));

% cube = zeros(183,183,90);
% cube(cst{4,4}{1}) = 1;
% mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin);
% 
% cst{12,1} = 11;
% cst{12,2} = 'Margin';
% cst{12,3} = 'OAR';
% cst{12,4}{1} = find(mVOIEnlarged);
% cst{12,5} = cst{11,5};
% 
% cst{12,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,30));
% 
% cube = zeros(183,183,90);
% cube(cst{10,4}{1}) = 1;
% mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin);
% 
% cst{13,1} = 12;
% cst{13,2} = 'Margin';
% cst{13,3} = 'OAR';
% cst{13,4}{1} = find(mVOIEnlarged);
% cst{13,5} = cst{11,5};
% 
% cst{13,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,30));
% 
% cube = zeros(183,183,90);
% cube(cst{6,4}{1}) = 1;
% cube(cst{1,4}{1}) = 1;
% mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin);
% 
% cst{14,1} = 13;
% cst{14,2} = 'Margin';
% cst{14,3} = 'OAR';
% cst{14,4}{1} = find(mVOIEnlarged);
% cst{14,5} = cst{11,5};
% 
% cst{14,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,30));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 120 240 270];
pln.propStf.couchAngles   = [0 0 0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% cube = zeros(183,183,90);
% sub2ind(pln.propStf.isoCenter(1,1),pln.propStf.isoCenter(1,2),pln.propStf.isoCenter(1,3));
% 
% dis = cst{4,4}{1} - pln.propStf.isoCenter(1,1);
% mir = pln.propStf.isoCenter(1,1) - dis;
% cube(mir) = 1;
% mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin);
% 
% cst{12,1} = 11;
% cst{12,2} = 'Margin';
% cst{12,3} = 'OAR';
% cst{12,4}{1} = find(mVOIEnlarged);
% cst{12,5} = cst{11,5};
% 
% cst{12,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,30));

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
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultGUImLETDose = resultGUI;

% cube = resultGUIRef.RBExD;
% plane = 3;
% slice = 34;
% doseWindow = [0 max(cube(:))];
% 
% figure,
% matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
% title('original plan')
% zoom(1.3)
% 
% cube = resultGUIRef.dirtyDose;
% plane = 3;
% slice = 34;
% doseWindow = [0 max(cube(:))];
% 
% figure,
% matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
% title('original plan')
% zoom(1.3)
% 
% maximumRef = max(resultGUIRef.RBExD(cst{9,4}{1}));
% minimumRef = min(resultGUIRef.RBExD(cst{9,4}{1}));

%% Proton Optimization 2 beam Body mean dirty Dose with penalty 100 dmax 20
clear("dij","pln","resultGUI","ct","cst","stf")

load("PROSTATE.mat")

cube = zeros(183,183,90);
cube(cst{6,4}{1}) = 1;
vResolution = ct.resolution;
vMargin = [];
vMargin.x = 5;
vMargin.y = 5;
vMargin.z = 5;
mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin);

cst{11,1}    = 10;
cst{11,2}    = 'Margin';
cst{11,3}    = 'OAR';
cst{11,4}{1} = find(mVOIEnlarged);
cst{11,5}    = cst{3,5};

cst{11,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(100,30)); 

cube = zeros(183,183,90);
cube(cst{4,4}{1}) = 1;
mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin);

cst{12,1} = 11;
cst{12,2} = 'Margin';
cst{12,3} = 'OAR';
cst{12,4}{1} = find(mVOIEnlarged);
cst{12,5} = cst{11,5};

cst{12,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,30));

cube = zeros(183,183,90);
cube(cst{10,4}{1}) = 1;
mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin);

cst{13,1} = 12;
cst{13,2} = 'Margin';
cst{13,3} = 'OAR';
cst{13,4}{1} = find(mVOIEnlarged);
cst{13,5} = cst{11,5};

cst{13,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,30));

cube = zeros(183,183,90);
cube(cst{6,4}{1}) = 1;
cube(cst{1,4}{1}) = 1;
mVOIEnlarged = matRad_addMargin(cube,cst,vResolution,vMargin);

cst{14,1} = 13;
cst{14,2} = 'Margin';
cst{14,3} = 'OAR';
cst{14,4}{1} = find(mVOIEnlarged);
cst{14,5} = cst{11,5};

cst{14,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(300,30));

cst{9,6}{2} = struct(DoseObjectives.matRad_SquaredOverdosing(300,maximumRef));
cst{6,6}{2} = struct(DirtyDoseObjectives.matRad_SquaredUnderdosingDirtyDose(10,10));

% Define radiation modality
pln.radiationMode = 'protons';        
pln.machine       = 'Generic';

% Calculate LET
pln.propDoseCalc.calcLET = 1;

% Set some beam parameters
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 120 240 270];
pln.propStf.couchAngles   = [0 0 0 0];
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
resultGUI= matRad_fluenceOptimization(dij,cst,pln);
resultTarget = resultGUI;

cube = resultTarget.RBExD;
plane = 3;
slice = 34;
% doseWindow = [0 max(cube(:))];

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
zoom(1.3)

cube = resultTarget.dirtyDose;
plane = 3;
slice = 34;
% doseWindow = [0 max(cube(:))];

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('manipulated plan')
zoom(1.3)

%% Differences

cube = resultGUIRef.RBExD - resultTarget.RBExD;
plane = 3;
slice = 34;
doseWindow = [-1 1];

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')
zoom(1.3)

cube = resultGUIRef.dirtyDose - resultTarget.dirtyDose;
plane = 3;
slice = 34;
% doseWindow = [0 max(cube(:))];

figure,
matRad_plotSliceWrapper(gca,ct,cst,1,cube,plane,slice,[],[],colorcube,[],doseWindow,[]);
title('difference plan')
zoom(1.3)
