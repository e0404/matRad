%% Example: Photon treatment plan with collapsed cone convolution dose calculation
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% In this example we will show
% (i) how to compile the collapsed cone transport mex file (optional)
% (ii) how to calculate the dose influence matrix of a photon plan with the
% collapsed cone convolution engine in a phantom with a lung insert
% (iii) how to recalculate the optimized plan with the SVD pencil beam
% engine and compare both dose distributions

%% Set matRad runtime configuration
matRad_rc;
matRad_cfg = MatRad_Config.instance();

%% Compile the transport (optional)
% Without the compiled mex file the engine uses a MATLAB implementation of
% the same transport, which gives the same result but is slower. With
% MATLAB on Linux and a compiler newer than the one supported by the MATLAB
% release, use matRad_compileCollapsedConeMex([], true, true) to link the
% C++ runtime statically.
if exist('matRad_collapsedConeMex', 'file') ~= 3
    try
        matRad_compileCollapsedConeMex();
    catch err
        matRad_cfg.dispWarning('Collapsed cone mex file not compiled: %s', err.message);
    end
end

%% Phantom with a lung insert
% A water phantom with a lung insert (relative electron density about 0.25)
% and a spherical tumor inside the lung.
builder = matRad_PhantomBuilder([80 80 50], [3 3 3], 1);
objTumor = struct(DoseObjectives.matRad_SquaredDeviation(800, 60));
objLung = struct(DoseObjectives.matRad_SquaredOverdosing(100, 20));
objBody = struct(DoseObjectives.matRad_SquaredOverdosing(10, 30));
builder.addSphericalTarget('Tumor', 6, 'objectives', objTumor, 'HU', 0);
builder.addBoxOAR('Lung', [36 40 30], 'objectives', objLung, 'HU', -750);
builder.addBoxOAR('Body', [76 76 46], 'objectives', objBody, 'HU', 0);
[ct, cst] = builder.getctcst();

%% Treatment plan
% The collapsed cone engine is selected with the short name 'CCC'. The
% number of cones, the scatter cut-off and the normalization depth can be
% set in pln.propDoseCalc.
pln.radiationMode = 'photons';
pln.machine = 'Generic';
pln.numOfFractions = 30;
pln.bioModel = 'none';
pln.multScen = 'nomScen';
pln.propStf.gantryAngles = 0:72:359;
pln.propStf.couchAngles = zeros(size(pln.propStf.gantryAngles));
pln.propStf.bixelWidth = 10;
pln.propStf.isoCenter = matRad_getIsoCenter(cst, ct, 0);
pln.propOpt.runDAO = false;
pln.propSeq.runSequencing = false;

pln.propDoseCalc.engine = 'CCC';
pln.propDoseCalc.numPolarAngles = 12;
pln.propDoseCalc.numAzimuthalAngles = 8;
pln.propDoseCalc.scatterCutOff = 150;
pln.propDoseCalc.doseGrid.resolution.x = 3;
pln.propDoseCalc.doseGrid.resolution.y = 3;
pln.propDoseCalc.doseGrid.resolution.z = 3;

%% Generate beam geometry
stf = matRad_generateStf(ct, cst, pln);

%% Dose influence matrix with the collapsed cone engine
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);

%% Fluence optimization
resultGUI = matRad_fluenceOptimization(dij, cst, pln);

%% Recalculation with the SVD pencil beam engine
plnSVD = pln;
plnSVD.propDoseCalc = struct('engine', 'SVDPB', 'doseGrid', pln.propDoseCalc.doseGrid);
resultSVD = matRad_calcDoseForward(ct, cst, stf, plnSVD, resultGUI.w);

%% Comparison
% The pencil beam algorithm scales the dose with the radiological depth
% only and overestimates the dose in and behind the lung for small fields,
% where the lateral electron equilibrium is not reached.
matRad_compareDose(resultGUI.physicalDose, resultSVD.physicalDose, ct, cst, [1 0 1], 'on', pln);

for i = 1:size(cst, 1)
    voxels = cst{i, 4}{1};
    matRad_cfg.dispInfo('%-6s mean dose CCC %6.2f Gy, SVDPB %6.2f Gy\n', cst{i, 2}, ...
                        mean(resultGUI.physicalDose(voxels)), mean(resultSVD.physicalDose(voxels)));
end
