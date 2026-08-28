%% Example: Neutron dose calculation (pencil beam kernels and MCNP export)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2025-2026 the matRad development team.
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
% (i)   how to set up a fast neutron therapy plan on a simple phantom with
%       different tissues (water, bone, lung),
% (ii)  how to calculate a neutron dose influence matrix with the SVD pencil
%       beam engine, which uses kernels fitted to Monte Carlo data and a
%       KERMA correction for the tissue types, and optimize the plan,
% (iii) how to export the plan as MCNP run files for an external Monte
%       Carlo calculation with the MCNP dose engine.

%% set matRad runtime configuration
matRad_rc; % If this throws an error, run it from the parent directory first to set the paths
matRad_cfg = MatRad_Config.instance();

%% Create a phantom with a target, a bone and a lung insert
% The phantom builder creates a water phantom; the inserts get their own HU
% values so that the KERMA correction of the neutron engine is visible.
ctDim = [60, 60, 40];        % x,y,z dimensions
ctResolution = [3, 3, 3];    % [mm], the MCNP engine needs equal x and y resolution

builder = matRad_PhantomBuilder(ctDim, ctResolution, 1);
builder.addBoxOAR('Body', ctDim .* ctResolution, 'HU', 0);
builder.addSphericalTarget('PTV', 15, 'offset', [0 -40 0], 'HU', 0, ...
                           'objectives', DoseObjectives.matRad_SquaredDeviation(800, 20));
builder.addSphericalOAR('Bone', 15, 'offset', [30 -10 0], 'HU', 1000, ...
                        'objectives', DoseObjectives.matRad_SquaredOverdosing(10, 10));
builder.addSphericalOAR('Lung', 15, 'offset', [-30 -10 0], 'HU', -700, ...
                        'objectives', DoseObjectives.matRad_SquaredOverdosing(10, 10));

[ct, cst] = builder.getctcst();

%% Treatment plan
% The 'Generic' neutron machine contains SVD kernels and a two-term depth
% dose parameterization fitted to MCNP simulations of the MEDAPP fast
% neutron therapy beam at FRM II, its neutron spectrum and HU-binned KERMA
% correction factors relative to water.
pln.radiationMode = 'neutrons';
pln.machine       = 'Generic';
pln.numOfFractions = 1;

pln.propStf.gantryAngles = [0 45 315];
pln.propStf.couchAngles  = [0 0 0];
pln.propStf.bixelWidth   = 10;

% Neutrons use the photon steering generators and the SVD pencil beam engine
pln.propDoseCalc.engine = 'SVDPB';
pln.propDoseCalc.doseGrid.resolution = struct('x', 3, 'y', 3, 'z', 3);

% Physical dose is optimized; a constant RBE can be used with
% pln.bioModel = 'constRBE' and pln.bioModel.RBE = 3 (Specht et al. 2015)
pln.bioModel = 'none';
pln.multScen = 'nomScen';

%% Generate beam geometry, calculate dose influence and optimize
stf = matRad_generateStf(ct, cst, pln);
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);
resultGUI = matRad_fluenceOptimization(dij, cst, pln);

%% Plot the result
% The lower dose in the bone insert (KERMA ratio to water of about 0.4)
% and the reduced attenuation in the lung insert are visible in the
% central slice.
slice = round(ct.cubeDim(3) / 2);
figure;
imagesc(resultGUI.physicalDose(:, :, slice));
axis equal tight; colorbar;
title('Neutron dose [Gy] in the central slice');

%% Export the plan for an MCNP calculation
% The MCNP engine writes one run file per bixel plus a script to run all of
% them into <userfolder>/MCNP/runfiles (externalCalculation = 'write').
% After running MCNP externally, set pln.propDoseCalc.externalCalculation
% to the folder containing the results to read them into a dij.
% Set exportForMCNP to true to write the run files (takes a while, as the
% CT is segmented into MCNP materials first).
exportForMCNP = false;

if exportForMCNP
    plnMCNP = pln;
    plnMCNP.propDoseCalc.engine = 'MCNP';
    plnMCNP.propDoseCalc.externalCalculation = 'write';
    plnMCNP.propDoseCalc.useDICOMinfoRescale = false;  % no DICOM rescale information for the phantom
    plnMCNP.propDoseCalc.bodyStructureName = 'Body';
    plnMCNP.propDoseCalc.lungStructureName = 'Lung';
    dijMCNP = matRad_calcDoseInfluence(ct, cst, stf, plnMCNP);
end
