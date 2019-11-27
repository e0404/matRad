% description
% 
% call
%    
%
% input
%   
%
% output 
%
%  
%
% References
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
load BOXPHANTOM_NARROW_NEW.mat

% meta information for treatment plan
pln.radiationMode   = 'protons';     % either photons / protons / carbon
pln.machine         = 'HITfixedBL';

pln.numOfFractions  = 30;
 
% beam geometry settings
pln.propStf.bixelWidth      = 50; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.longitudinalSpotSpacing = 50;
pln.propStf.gantryAngles    = 0; % [?] 
pln.propStf.couchAngles     = 0; % [?]
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = [ct.cubeDim(2) / 2 * ct.resolution.y, ...
                                0, ...
                                ct.cubeDim(3) / 2 * ct.resolution.z];

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

%% baseData fitting
 
%read MC phase space data
dataMC = readtable('BDL_matrad.txt');             
dataMC = dataMC{11:end,:};
tmp = [];
for i = 1:size(dataMC,1)    
    tmp = [tmp; strsplit(dataMC{i})];
end
energyMC    = str2double(tmp(:,1));
spotMC      = str2double(tmp(:,6));
divMC       = str2double(tmp(:,7));
corMC       = str2double(tmp(:,8));
  
minEnergy = 70;
maxEnergy = 225;
nEnergy   = 100;  

count = 1;
tic
for currentEnergy = linspace(minEnergy, maxEnergy, nEnergy)
    
    % calculate sigma/spotsize at isocenter using MC phase space data
    divNozzle = interp1(energyMC,divMC,currentEnergy); 
    corNozzle = interp1(energyMC,corMC,currentEnergy); 
    spotNozzle = interp1(energyMC,spotMC,currentEnergy); 
    z = 420;
    spotIso = sqrt(spotNozzle^2 + 2 * (corNozzle*spotNozzle + divNozzle * z) * divNozzle * z - divNozzle^2*z^2);
    
    %assign energy to stf and run MC simulation
    stf.ray.energy = currentEnergy;
    resultGUI = matRad_calcDoseDirectMC(ct,stf,pln,cst,ones(sum(stf(:).totalNumOfBixels),1),10000000);          
    
    machine.data(count) = matRad_fitBaseData(resultGUI.physicalDose, ct.resolution, currentEnergy, spotIso);

    disp(['baseData Progress :', ' ', num2str(round(count/nEnergy*100)), '%']);
    count = count + 1;
end
toc

%save data in machine
machine.meta.radiationMode = 'protons';
machine.meta.dataType = 'singleGauss';
machine.meta.created_on = date;
machine.meta.created_by = 'Paul Anton Meder';
machine.meta.SAD = (2218 + 1839) / 2;
machine.meta.BAMStoIsoDist = 420.0;
machine.meta.machine = 'Generic';
machine.meta.LUT_bxWidthminFWHM = [1, Inf; 8 ,8];
  

datetime