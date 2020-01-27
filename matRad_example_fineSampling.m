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

%load HEAD_AND_NECK
% load TG119.mat
% load slabPhantom6.mat
% load LungPatient1.mat
% load PROSTATE.mat
% load LIVER.mat
load BOXPHANTOM
% load BOXPHANTOMv3.mat
% load BOXPHANTOM_NARROW_NEW.mat
% load phantomTest.mat


% meta information for treatment plan
pln.radiationMode   = 'protons';     % either photons / protons / carbon
pln.machine         = 'matRadBDLfixed';

pln.numOfFractions  = 30;

% beam geometry settings
pln.propStf.bixelWidth      = 300; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.longitudinalSpotSpacing = 300;
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
load protons_matRadBDL;
stf.ray.energy = machine.data(44).energy;
                            
%% dose calculation
 % analytical dose without fine sampling
    tic
    pln.propDoseCalc.anaMode = 'standard';
    dij = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI = matRad_calcCubes(ones(sum(stf(:).totalNumOfBixels),1),dij);
    anaDose     = resultGUI.physicalDose;
    t1 = toc;

%  % analytical dose with fine sampling
    tic
    pln.propDoseCalc.anaMode = 'fineSampling';
    dijFS = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI_FS = matRad_calcCubes(ones(sum(stf(:).totalNumOfBixels),1),dijFS);
    resultGUI.physicalDoseFS = resultGUI_FS.physicalDose;
    anaFsDose   = resultGUI.physicalDoseFS;
    t2 = toc;

 % analytical dose with std correction
    tic
    pln.propDoseCalc.anaMode = 'stdCorr';
    dijSC = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI_SC = matRad_calcCubes(ones(sum(stf(:).totalNumOfBixels),1),dijSC);
    resultGUI.physicalDoseSC = resultGUI_SC.physicalDose;
    anaScDose   = resultGUI.physicalDoseSC;
    t3 = toc;

 % Monte Carlo dose
    tic
    resultGUI_MC = matRad_calcDoseDirectMC(ct,stf,pln,cst,ones(sum(stf(:).totalNumOfBixels),1), 10000000);
    resultGUI.physicalDoseMC = resultGUI_MC.physicalDose;
    mcDose      = resultGUI.physicalDoseMC;
    t4 = toc;

 %% plot doses
figure
subplot(2,2,1)
imagesc(anaDose(38:122,38:122,round(ct.cubeDim(3)/2)));
caxis([0 2e-3]);
title('Standard dose');
hold on 
contour(ct.cube{1}(38:122,38:122,round(ct.cubeDim(3)/2)),3,'color','white');
hold off

subplot(2,2,2)
imagesc(anaFsDose(38:122,38:122,round(ct.cubeDim(3)/2)));
caxis([0 2e-3]);
title('Fine sampling dose');
hold on
contour(ct.cube{1}(38:122,38:122,round(ct.cubeDim(3)/2)),3,'color','white');
hold off

subplot(2,2,3)
imagesc(mcDose(38:122,38:122,round(ct.cubeDim(3)/2)));
caxis([0 2e-3]);
title('Monte Carlo dose');
hold on 
contour(ct.cube{1}(38:122,38:122,round(ct.cubeDim(3)/2)),3,'color','white');
hold off

subplot(2,2,4)
imagesc(anaScDose(38:122,38:122,round(ct.cubeDim(3)/2)));
caxis([0 2e-3]);
title('Std corrected dose');
hold on 
contour(ct.cube{1}(38:122,38:122,round(ct.cubeDim(3)/2)),3,'color','white');
hold off

 %% execute gamma tests
gammaTest = [2, 2];
figure
subplot(2,2,3.5)
[gammaCube1,gammaPassRateCell] = matRad_gammaIndex(anaDose,mcDose,[ct.resolution.x,ct.resolution.y,ct.resolution.z],gammaTest,round(ct.cubeDim(3)/2),3,'global',cst);
title({[num2str(gammaPassRateCell{1,2}) '% of points > ' num2str(gammaTest(1)) '% pass gamma criterion (' num2str(gammaTest(1)) '%/ ' num2str(gammaTest(2)) 'mm)'];'stadard dose'});
% axis([38 122 38 122]);
hold on
contour(ct.cube{1}(:,:,round(ct.cubeDim(3)/2)),3,'color','white');
hold off

subplot(2,2,2)
[gammaCube2,gammaPassRateCell] = matRad_gammaIndex(anaFsDose,mcDose,[ct.resolution.x,ct.resolution.y,ct.resolution.z],gammaTest,round(ct.cubeDim(3)/2),3,'global',cst);
title({[num2str(gammaPassRateCell{1,2}) '% of points > ' num2str(gammaTest(1)) '% pass gamma criterion (' num2str(gammaTest(1)) '%/ ' num2str(gammaTest(2)) 'mm)'];'fine sampling dose'});
% axis([38 122 38 122]);
hold on
contour(ct.cube{1}(:,:,round(ct.cubeDim(3)/2)),3,'color','white');
hold off

subplot(2,2,1)
[gammaCube3,gammaPassRateCell] = matRad_gammaIndex(anaScDose,mcDose,[ct.resolution.x,ct.resolution.y,ct.resolution.z],gammaTest,round(ct.cubeDim(3)/2),3,'global',cst);
title({[num2str(gammaPassRateCell{1,2}) '% of points > ' num2str(gammaTest(1)) '% pass gamma criterion (' num2str(gammaTest(1)) '%/ ' num2str(gammaTest(2)) 'mm)'];'std corrected dose'});
% axis([38 122 38 122]);
hold on
contour(ct.cube{1}(:,:,round(ct.cubeDim(3)/2)),3,'color','white');
hold off
