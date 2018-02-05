%% prep

%clearvars -except *dir
%close all

% load patient data, i.e. ct, voi, cst

%load HEAD_AND_NECK
load TG119.mat
%load PROSTATE.mat
%load LIVER.mat
%load BOXPHANTOM.mat

% meta information for treatment plan
pln.isoCenter       = matRad_getIsoCenter(cst,ct,0);
pln.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.gantryAngles    = [0:72:359]; % [°]
pln.couchAngles     = [0 0 0 0 0]; % [°]
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = prod(ct.cubeDim);
pln.voxelDimensions = ct.cubeDim;
pln.radiationMode   = 'photons';     % either photons / protons / carbon
pln.bioOptimization = 'none';        % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                     % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
%pln.numOfFractions  = 38;
pln.runSequencing   = true; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.runDAO          = true; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.VMAT            = false; % 1/true: run VMAT, 0/false: don't
pln.dynamic         = true;
pln.halfFluOpt      = false; % indicates if you want to constrain half of each field to have 0 fluence (other half is compensated on the other side)
pln.machine         = 'Generic';

% For VMAT
pln.runSequencing   = true;
pln.runDAO          = true;
pln.VMAT            = true;

pln.scaleDRx        = true;
pln.scaleDij        = true;
pln.jacobi          = true;
pln.dynamic         = true;


pln.numApertures = 7; %max val is pln.maxApertureAngleSpread/pln.minGantryAngleRes
pln.numLevels = 7;

pln.minGantryAngleRes = 4; %Bzdusek
pln.maxApertureAngleSpread = 28; %should be an even multiple of pln.minGantryAngleRes; Bzdusek

pln = matRad_VMATGantryAngles(pln,'new');


pln.gantryRotCst = [0 6]; %degrees per second
pln.defaultGantryRot = max(pln.gantryRotCst); %degrees per second
pln.leafSpeedCst = [0 6]*10; %mm per second
pln.defaultLeafSpeed = pln.leafSpeedCst(2);
pln.doseRateCst = [75 600]/60; %MU per second
pln.defaultDoseRate = pln.doseRateCst(2);
pln.maxLeafTravelPerDeg = pln.leafSpeedCst(2)/pln.defaultGantryRot;

pln.halfFluOpt = false;
pln.halfFluOptMargin = 10; % mm

recalc.doRecalc = 0;
recalc.dynamic = 0;
recalc.interpNew = 0;
recalc.pln = pln;
recalc.pln.minGantryAngleRes = 1;


if ~exist('stf','var')
    stf = matRad_generateStf(ct,cst,pln);
end

%load Dij

if ~exist('dij','var')
    dij = loadDij('TG119');
end

% inverse planning for imrt
resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf,0);

%% run DAO 4 times with different sequencing algorithms

% Siochi
fname = 'Siochi50';
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);

t0_nDij_nJ = tic;
resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln);
t_nDij_nJ = toc(t0_nDij_nJ);
savefig(fname)
save(fname,'resultGUI')

%{
% Svensson
fname = 'Svensson';
resultGUI = matRad_svenssonLeafSequencing(resultGUI,stf,dij,pln,0);

t0_yDij_nJ = tic;
resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln);
t_yDij_nJ = toc(t0_yDij_nJ);
savefig(fname)
save(fname,'resultGUI')
%}

save('timings','t_*')




