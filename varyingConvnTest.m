% matRad script
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_rc

% load patient data, i.e. ct, voi, cst

% load PHANTOM_control.mat
% load PHANTOM_hetSlab_entrance_10mm.mat
% load PHANTOM_hetSlab_entrance_30mm.mat
load PHANTOM_slab_entrance_10mm.mat
% load PHANTOM_slab_entrance_50mm.mat
% load PHANTOM_slab_entrance_100mm.mat
% load PHANTOM_slab_mid_10mm.mat
% load PHANTOM_wedge_entrance_30mm.mat
% load PHANTOM_wedge_entrance_100mm.mat
% load PHANTOM_wedge_mid_30mm.mat

% load LungPatient1.mat

% meta information for treatment plan
pln.radiationMode   = 'protons';     % either photons / protons / carbon
pln.machine         = 'matRadBDL_APM';

pln.numOfFractions  = 30;

% beam geometry settings
pln.propStf.bixelWidth      = 500; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.longitudinalSpotSpacing = 500;
pln.propStf.gantryAngles    = 0; % [?] 
pln.propStf.couchAngles     = 0; % [?]
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
                            
% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = ct.resolution.x; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = ct.resolution.y; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = ct.resolution.z; % [mm]
pln.propDoseCalc.airOffsetCorrection = true;

% optimization settings
pln.propOpt.optimizer       = 'IPOPT';
pln.propOpt.bioOptimization = 'none'; % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                      % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;  % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

%% generate steering file
stf = matRad_generateStf(ct,cst,pln);

% assign new energy
load protons_matRadBDL_APM;
stf.ray.energy = machine.data(36).energy;

%% dose calculation
% analytical standard dose
pln.propDoseCalc.anaMode = 'standard';
dij = matRad_calcParticleDose(ct,stf,pln,cst,false);
resultGUI = matRad_calcCubes(ones(sum(stf(:).totalNumOfBixels),1),dij);
anaDose     = resultGUI.physicalDose;
    
% analytical dose with std correction
pln.propDoseCalc.anaMode = 'stdCorr';
dijSC = matRad_calcParticleDose(ct,stf,pln,cst,false,1);
resultGUI_SC = matRad_calcCubes(ones(sum(stf(:).totalNumOfBixels),1),dijSC);
resultGUI.physicalDoseSC = resultGUI_SC.physicalDose;
anaScDose   = resultGUI.physicalDoseSC;

%% plotting
contourSwitch = true;
figure;
subplot(1,2,1)
imagesc(anaDose(:,:,round(ct.cubeDim(3)/2)));
caxis([0 2e-3]);
title('Standard dose');
hold on 
contour(ct.cube{1}(:,:,round(ct.cubeDim(3)/2)),1,'color','white');
if contourSwitch
    contour(anaDose(:,:,round(ct.cubeDim(3)/2)),10,'color','black');
end
hold off
pbaspect([ct.cubeDim(2) ct.cubeDim(1) ct.cubeDim(3)])

subplot(1,2,2)
imagesc(anaScDose(:,:,round(ct.cubeDim(3)/2)));
caxis([0 2e-3]);
title('Std corrected dose');
hold on 
contour(ct.cube{1}(:,:,round(ct.cubeDim(3)/2)),1,'color','white');
if contourSwitch
    contour(anaScDose(:,:,round(ct.cubeDim(3)/2)),10,'color','black');
end
hold off
pbaspect([ct.cubeDim(2) ct.cubeDim(1) ct.cubeDim(3)])
