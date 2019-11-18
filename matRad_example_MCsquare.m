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
pln.machine         = 'HITmcSquare';

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

load protons_HITmcSquare.mat
stf.ray.energy = machine.data(1).energy;
% 
% %% meanEnergy vs nomEnergy 
% p1 =  -5.472e-05;
% p2 =        1.03;
% p3 =      -3.109;
% meanEnergy = @(x) p1*x^2 + p2*x + p3;
% 
% %% spreads vs nomEnergy
% q1 =  -6.406e-12;
% q2 =   4.764e-09;
% q3 =  -1.452e-06;
% q4 =   0.0002221;
% q5 =    -0.01699;
% q6 =      0.9305;
%  
% spread = @(x) q1*x^5 + q2*x^4 + q3*x^3 + q4*x^2 + q5*x + q6;
      
       
%% dose calculation
if strcmp(pln.radiationMode,'photons')
    dij = matRad_calcPhotonDose(ct,stf,pln,cst);
    dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
    dij = matRad_calcParticleDose(ct,stf,pln,cst);
    dijMC = matRad_calcParticleDoseMC(ct,stf,pln,cst,1000000);
end

resultGUI = matRad_calcCubes(ones(dij.totalNumOfBixels,1),dij);
resultGUI_MC = matRad_calcCubes(resultGUI.w,dijMC);

resultGUI.physicalDose_MC = resultGUI_MC.physicalDose;
mcIDD  = sum(sum(resultGUI.physicalDose_MC,2),3);
anaIDD = sum(sum(resultGUI.physicalDose,2),3);
sqDev = sum((mcIDD - anaIDD).^2);

plot(mcIDD)
hold on
plot(anaIDD)
resultGUI.physicalDose_diff = (resultGUI.physicalDose - resultGUI.physicalDose_MC);

matRadGUI;



% 
% %% write new machine
% 
% for i = 1:numel(machine.data)
%     for k = 1:4
%         machine.data(i).mcSquareData(k).MeanEnergy   = meanEnergy(machine.data(i).mcSquareData(k).NominalEnergy);
%         machine.data(i).mcSquareData(k).EnergySpread = spread(machine.data(i).mcSquareData(k).NominalEnergy);
%     end
% end









%%                       
% num = 22;
% spreads       = linspace( 0, 0.6, num);
% energyOffsets = linspace(-3, 2.5, num);
% 
% count = 1;
% for centerEnergy = [machine.data(1:18:end).energy, machine.data(end).energy]
%     
%     erg{count,1} = centerEnergy;
%     
%     stf.ray.energy = centerEnergy;
%     dij = matRad_calcParticleDose(ct,stf,pln,cst);
%     resultGUI = matRad_calcCubes(ones(dij.totalNumOfBixels,1),dij);
%     ixE = 1;
%     for offset = energyOffsets
%     
%         mean = centerEnergy + offset;
%         erg{count, 2}(1,ixE) = mean;
%         
%         ixSpread = 1;
%         for spread = spreads
%             erg{count, 3}(1,ixSpread) = spread;
%             dijMC = matRad_calcParticleDoseMC(ct,stf,pln,cst,1000000,0,mean,spread);
%             resultGUI_MC = matRad_calcCubes(resultGUI.w,dijMC);
%             resultGUI.physicalDose_MC = resultGUI_MC.physicalDose;
% 
%             mcIDD  = sum(sum(resultGUI.physicalDose_MC,2),3);
%             anaIDD = sum(sum(resultGUI.physicalDose,2),3);
%             erg{count,4}(ixE,ixSpread) = sum((mcIDD - anaIDD).^2);
%             
%             
%             ixSpread = ixSpread + 1;
%         end       
%         ixE = ixE + 1;
%     end
%     count = count + 1;
% end

