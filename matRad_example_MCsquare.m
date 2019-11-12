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
pln.machine         = 'generic_MCsquare';

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
load protons_generic_MCsquare
stf.ray.energy = machine.data(end).energy;

%% dose calculation
if strcmp(pln.radiationMode,'photons')
    dij = matRad_calcPhotonDose(ct,stf,pln,cst);
    %dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
    dij = matRad_calcParticleDose(ct,stf,pln,cst);
    dijMC = matRad_calcParticleDoseMC(ct,stf,pln,cst,1000000);
end

resultGUI = matRad_calcCubes(ones(dij.totalNumOfBixels,1),dij);
resultGUI_MC = matRad_calcCubes(resultGUI.w,dijMC);

resultGUI.physicalDose_MC = resultGUI_MC.physicalDose;
resultGUI.physicalDose_diff = (resultGUI.physicalDose - resultGUI.physicalDose_MC);

matRadGUI;




                      
% num = 100;
% array = linspace(0,4,num);
% 
% count = 1;
% for e = [machine.data(1:4:end).energy]
% 
%     erg{count, 1} = e;
%     stf.ray.energy = e;
%     [~ ,i, ~] = intersect([machine.data(:).energy], e);
%     
%     dij = matRad_calcParticleDose(ct,stf,pln,cst);
%     resultGUI = matRad_calcCubes(ones(dij.totalNumOfBixels,1),dij);
% 
%     newDepths = linspace(0,machine.data(i).depths(end),numel(machine.data(i).depths) * 100);
%     newDose   = interp1(machine.data(i).depths, machine.data(i).Z, newDepths, 'spline');
%      
%     machine.data(i).depths = newDepths;
%     machine.data(i).Z      = newDose;
%     
%     %interpolate range at 80% dose after peak.
%     [maxV, maxI] = max(machine.data(i).Z);
%     [~, r80ind] = min(abs(machine.data(i).Z(maxI:end) - 0.8 * maxV));
%  
%     %find FWHM w50 of bragg peak
%     [~, d50rInd] = min(abs(machine.data(i).Z(maxI:end) - 0.5 * maxV));
%     d50rInd = d50rInd - 1;
%     d50_r = interp1(machine.data(i).Z(maxI + d50rInd - 1:maxI + d50rInd + 1), ...
%                             machine.data(i).depths(maxI + d50rInd - 1:maxI + d50rInd + 1), 0.5 * maxV);
%                              
%     array1 = machine.data(i).depths(maxI-round(d50rInd):maxI+d50rInd);
%     array2 = machine.data(i).Z(maxI-round(d50rInd):maxI+d50rInd);  
% 
%     gauss = @(x, a, sigma, b) a * exp(- (x - b).^2 / (2*sigma^2));
% 
%     funcs.objective = @(p) sum((gauss(array1, p(1), p(2),p(3)) - array2).^2);
% 
%     funcs.gradient = @(p) [ sum(2 * (p(1) * exp(- (array1 - p(3)).^2 / (2 * p(2)^2)) - array2) ...
%                                             .* exp(- (array1 - p(3)).^2 / (2 * p(2)^2)));
%                             sum(2 * (p(1) * exp(- (array1 - p(3)).^2 / (2 * p(2)^2)) - array2) ...
%                                             .* p(1) .* exp(- (array1 - p(3)).^2 / (2 * p(2)^2)) .* (array1 - p(3)).^2 / p(2)^3)
%                             sum(2 * (p(1) * exp(- (array1 - p(3)).^2 / (2 * p(2)^2)) - array2) ...
%                                             .* p(1) .* exp(- (array1 - p(3)).^2 / (2 * p(2)^2)) .* 2 .* (array1 - p(3)) / (2 * p(2)^2))];
% 
% 
%     options.lb = [0,  0, 0];
%     options.ub = [ Inf, Inf, Inf];
%     options.ipopt.limited_memory_update_type = 'bfgs';
%             options.ipopt.hessian_approximation = 'limited-memory';
%     options.ipopt.print_level = 1;
% 
%     start = [maxV, 2 * d50_r, machine.data(i).depths(maxI)];
%     [fitResult, ~] = ipopt (start, funcs, options);
%     
%     
%     erg{count, 2} = fitResult(2);
%     deviation = [];
%     
%     for spread = array
% 
%         dijMC = matRad_calcParticleDoseMC(ct,stf,pln,cst,1000000,0,spread);
%         resultGUI_MC = matRad_calcCubes(resultGUI.w,dijMC);
%         resultGUI.physicalDose_MC = resultGUI_MC.physicalDose;
%         resultGUI.physicalDose_diff = (resultGUI.physicalDose - resultGUI.physicalDose_MC);
%         
%         deviation = [deviation, sum(sum(sum(resultGUI.physicalDose_diff.^2)))];
%     end
%     
%     erg{count, 3} = deviation;
%     count = count + 1;
% end
% 
% %%
% num = 100;
% array = linspace(0,4,num);
% 
%      
% 
% tmp = [];
% for i = 1:19
% %     plot(erg{i,3})
%     [~, index] = min(erg{i,3})
%     tmp = [tmp, array(index)];
% end
% 
% scatter([erg{1:19,2}], tmp)
% ylim([0, 0.5]);

