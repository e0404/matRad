%% Example: Proton Treatment Plan with FRED MC
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% In this example we will show 
% (i)  how to compute a simple plan using the FRED MC engine
% (ii) how to compute LETd distributions and apply a biological model

%% Setup the plan parameters
matRad_rc;
matRad_cfg = MatRad_Config.instance();

load('TG119.mat');

pln.radiationMode = 'protons';
pln.machine       = 'Generic';

pln.propDoseCalc.calcLET = 0;

pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [30,330];
pln.propStf.couchAngles   = zeros(size(pln.propStf.gantryAngles));
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propSeq.runSequencing = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

% Start with a simple physical dose model
pln.bioParam = matRad_bioModel(pln.radiationMode,'none');

pln.multScen = matRad_multScen(ct, 'nomScen');

stf = matRad_generateStf(ct,cst,pln);

%% Let's start with the analytical dose calculation algorithm
pln.propDoseCalc.engine = 'HongPB';

% Compute the dose influence matrix
dij_PB = matRad_calcDoseInfluence(ct,cst,stf,pln);

%% Fluence optimization
resultGUI_PB = matRad_fluenceOptimization(dij_PB,cst,pln);
wOptimized = resultGUI_PB.w;

%% Let's now re-compute the dose distribution with FRED MC
pln.propDoseCalc.engine = 'FRED';
pln.propDoseCalc.useGPU = false;

% For illustraton, let's compute a dose influence matrix.
dij_FRED = matRad_calcDoseInfluence(ct,cst,stf,pln);

% And compute the dose cube from the dij:
resultGUI_FRED = matRad_calcCubes(wOptimized,dij_FRED,1);

% Alternative is to perform a direct dose calculation:
% resultGUI_FRED = matRad_calcDoseForward(ct,cst,stf,pln,wOptimized)

%% Compare doses
matRad_compareDose(resultGUI_PB.physicalDose, resultGUI_FRED.physicalDose,ct,cst,[],'on');
%% Tune the MC calculation parameters
% This example illustrates some of the parameters that can be tuned for the
% Engine. For more information run:
% help DoseEngines.matRad_ParticleFREDEngine

pln.bioParam = matRad_bioModel(pln.radiationMode,'MCN');

pln.propDoseCalc.useGPU      = true;
pln.propDoseCalc.sourceModel = 'gaussian'; %alternatives: {'gaussian', 'emittance', 'sigmaSqrModel'}
pln.propDoseCalc.HUtable     = 'internal'; % default: matRad_default_FredMaterialConverter
pln.propDoseCalc.scorers     = {'Dose', 'LETd'};

resultGUI_recalc = matRad_calcDoseForward(ct,cst,stf,pln, wOptimized);

%% Compare physical dose and RBExD distributions
matRad_compareDose(resultGUI_recalc.physicalDose, resultGUI_recalc.RBExDose,ct,cst,[],'on');