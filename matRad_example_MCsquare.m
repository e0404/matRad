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
%load TG119.mat
%load PROSTATE.mat
%load LIVER.mat
load BOXPHANTOMv3.mat

% meta information for treatment plan
pln.radiationMode   = 'protons';     % either photons / protons / carbon
pln.machine         = 'HITfixedBL';

pln.numOfFractions  = 30;

% beam geometry settings
pln.propStf.bixelWidth      = 50; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.longitudinalSpotSpacing = 100;
pln.propStf.gantryAngles    = 0; % [?] 
pln.propStf.couchAngles     = 0; % [?]
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = ct.resolution.x; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = ct.resolution.y; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = ct.resolution.z; % [mm]

% optimization settings
pln.propOpt.optimizer       = 'IPOPT';
pln.propOpt.bioOptimization = 'none'; % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                      % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;  % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

%% generate steering file
stf = matRad_generateStf(ct,cst,pln);
load protons_HITfixedBL
stf.ray.energy = machine.data(120).energy;

eSpread = @(E)0.3566*(1-exp(-0.03307*(E-57.25))) .* heaviside(E-57.25);



%% dose calculation
if strcmp(pln.radiationMode,'photons')
    dij = matRad_calcPhotonDose(ct,stf,pln,cst);
    %dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
    dij = matRad_calcParticleDose(ct,stf,pln,cst);
    dijMC = matRad_calcParticleDoseMC(ct,stf,pln,cst,100000000,0,eSpread(stf.ray.energy));
end

resultGUI = matRad_calcCubes(ones(dij.totalNumOfBixels,1),dij);
resultGUI_MC = matRad_calcCubes(resultGUI.w,dijMC);


resultGUI.physicalDose_MC = resultGUI_MC.physicalDose;
resultGUI.physicalDose_diff = (resultGUI.physicalDose - resultGUI.physicalDose_MC);


matRadGUI;

% num = 102;
% 
% array = linspace(0,0.9,num);
% sqDev = zeros(numel(machine.data), num);
% tic
% i = 1;
% for e = [machine.data(1:6:end).energy]
% 
%     stf.ray.energy = e;
%     
%     dij = matRad_calcParticleDose(ct,stf,pln,cst);
%     resultGUI = matRad_calcCubes(ones(dij.totalNumOfBixels,1),dij);
% 
%     count = 1;
%     
%     for spread = array
%         i ,count
% 
%         dijMC = matRad_calcParticleDoseMC(ct,stf,pln,cst,1000000,0,spread);
% 
%         resultGUI_MC = matRad_calcCubes(resultGUI.w,dijMC);
% 
%         resultGUI.physicalDose_MC = resultGUI_MC.physicalDose;
% 
%         resultGUI.physicalDose_diff = (resultGUI.physicalDose - resultGUI.physicalDose_MC);
% 
%         sqDev(i, count) = sum(sum(sum(resultGUI.physicalDose_diff.^2)));
%         count = count + 1;
%     end
%     
%     i = i + 1;
% end
% toc
% 
%  %%
% best = [];
% for t = 1:size(actual, 1)
%     
%     % % plot(array, sqDev);
%     p = polyfit(array, actual(t,:), 2);
%     %p = @(x) p(1)*x.^2 + p(2)*x + p(3);
% 
% %     plot(array, actual(t,:))
% %     hold on
% %     plot(array, p(array))
% %     hold off
%     best = [best, p(2)/(-2*p(1))];
% 
%     
% end
% % matRadGUI;
