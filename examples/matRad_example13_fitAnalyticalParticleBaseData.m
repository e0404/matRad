%% Example: baseData fitting to mcSquare simulation
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %% 1) CT & CST creation
clear
matRad_rc

ixOAR = 1;
ixPTV = 2;

% define general VOI properties
cst{ixOAR,1} = 0;
cst{ixOAR,2} = 'contour';
cst{ixOAR,3} = 'OAR';

cst{ixPTV,1} = 1;
cst{ixPTV,2} = 'target';
cst{ixPTV,3} = 'TARGET';

% define optimization parameter for both VOIs
cst{ixOAR,5}.TissueClass = 1;
cst{ixOAR,5}.alphaX      = 0.1000;
cst{ixOAR,5}.betaX       = 0.0500;
cst{ixOAR,5}.Priority    = 2;
cst{ixOAR,5}.Visible     = 1;
cst{ixOAR,5}.visibleColor     = [1 0 0];
cst{ixOAR,6}{1,1}.className   = 'DoseObjectives.matRad_SquaredOverdosing';
cst{ixOAR,6}{1,1}.parameters{1}  = 5;
cst{ixOAR,6}{1,1}.penalty     = 100;


cst{ixPTV,5}.TissueClass = 1;
cst{ixPTV,5}.alphaX      = 0.1000;
cst{ixPTV,5}.betaX       = 0.0500;
cst{ixPTV,5}.Priority    = 1;
cst{ixPTV,5}.Visible     = 1;
cst{ixPTV,5}.visibleColor     = [0 1 0];
cst{ixPTV,6}{1,1}.className   = 'DoseObjectives.matRad_SquaredOverdosing';
cst{ixPTV,6}{1,1}.parameters{1}  = 60;
cst{ixPTV,6}{1,1}.penalty     = 800;


%% Create CT
xDim = 1100;
yDim = 250;
zDim = 250;

cubeDim      = [xDim yDim zDim];
ct.cubeDim = cubeDim;

ct.cube{1} = ones(cubeDim) * 1;
ct.cube{1}(1,1,1) = 0; 

ct.resolution.x = 0.32;
ct.resolution.y = 0.32;
ct.resolution.z = 0.32;

ct.numOfCtScen  = 1;


%% Create a cubic phantom
iso = [800,125,125];

% create an ct image series with zeros
ct.cubeHU{1} = ones(ct.cubeDim) * 0;
ct.cubeHU{1}(1,1,1) = -1000; 

ct.hlut = [1,0;0,-1024];

% create body of full phantom size
body = ones(ct.cubeDim);
cst{1,4}{1} = find(body == 1);

% create target
centerP_corr = iso;
height_corr = 10;
width_corr = 10;
depth_corr = 10;

target = zeros(ct.cubeDim);
for i=-height_corr/2+1:height_corr/2
    for j=-width_corr/2:width_corr/2
        for k=-depth_corr/2:depth_corr/2
            target(centerP_corr(1)+i,centerP_corr(2)+j,centerP_corr(3)+k) = 1;
        end
    end
end
cst{2,4}{1} = find(target == 1);

disp('CT creation done!');
clearvars -except ct cst matRad_cfg


%% 2) MCsquare computation and baseData fitting

% meta information for treatment plan
pln.radiationMode   = 'protons'; 

% create meta machine data
machine.meta.machine = 'example'; %name of the machine
machine.meta.radiationMode = 'protons'; %modality
machine.meta.dataType = 'singleGauss'; %singleGauss or doubleGauss
machine.meta.created_on = date;
machine.meta.created_by = 'matRad_example';
machine.meta.SAD = (2218 + 1839) / 2; %This is the (virtual) source to axis distance
machine.meta.BAMStoIsoDist = 420.0; %distance from beam nozzle ot isocenter
machine.meta.LUT_bxWidthminFWHM = [0, Inf; 5 ,5]; %Specifies which minimum FWHM to use as spot sice for which ranges of lateral spot distance (here, each spot distance of 0 to to Inf gets at least 5mm wide spots
machine.meta.fitAirOffset = 420.0; %Tells matRad how much "air" was considered during fitting. Set this to 0 if the fit is obtained in vacuum and no air transport is simulated up to the phantom. matRad assumes that the phantom starts at the isocenter.

% Now add the example machine to the pln and then save it
pln.machine = machine.meta.machine;
pln.radiationMode = machine.meta.radiationMode;
fileName = [pln.radiationMode '_' pln.machine];
filePath = fullfile(matRad_cfg.matRadRoot,'basedata',[fileName '.mat']);

matRad_cfg.dispInfo('Saving temporary machine %s to %s\n',fileName,filePath);
save(filePath,'machine','-v7');
clear machine; 

 
% beam geometry settings
pln.propStf.bixelWidth      = 50; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.longitudinalSpotSpacing = 50;
pln.propStf.gantryAngles    = 0; % [?] 
pln.propStf.couchAngles     = 0; % [?]
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);

% set isoCenter at entrance of beam into phantom
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

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen');


quantityOpt   = 'physicalDose';            % either  physicalDose / effect / RBExD
modelName     = 'none';         % none: for photons, protons, carbon                                    constRBE: constant RBE model
                                    % MCN: McNamara-variable RBE model for protons                          WED: Wedenberg-variable RBE model for protons 
                                    % LEM: Local Effect Model for carbon ions
% retrieve bio model parameters
pln.bioParam = matRad_BioModel(pln.radiationMode,quantityOpt, modelName);


%% generate steering file
stf = matRad_generateStfSinglePencilBeam(ct,cst,pln);

% Select existing BDL file to load and fit
pln.loadExistingBDL = 'BDL_matRad.txt';

% read MC phase space data
dataMC = importdata(pln.loadExistingBDL, ' ', 16);       
energyMC = dataMC.data(:, 1);
spotMC   = (dataMC.data(:, 6) + dataMC.data(:, 9)) / 2;
divMC    = (dataMC.data(:, 7) + dataMC.data(:,10)) / 2;
corMC    = (dataMC.data(:, 8) + dataMC.data(:,11)) / 2;

%% Run Base Data Fitting

% define energy range
minEnergy = 70;
maxEnergy = 225;
nEnergy   = 75;

% Number of histories for the MC simulation
pln.propMC.numHistories = 1e5;

% We create a figure to display the fit
hf = figure();

% Here we loop over all energies we want to fit
count = 1;
for currentEnergy = linspace(minEnergy, maxEnergy, nEnergy)
    
    % calculate sigma/spotsize at isocenter using MC phase space data
    divNozzle = interp1(energyMC,divMC,currentEnergy); 
    corNozzle = interp1(energyMC,corMC,currentEnergy); 
    spotNozzle = interp1(energyMC,spotMC,currentEnergy); 
    z = 420;
    
    mcData.divNozzle = divNozzle;
    mcData.corNozzle = corNozzle;
    mcData.spotNozzle = spotNozzle;
    mcData.z = z;
    spotIso = sqrt(spotNozzle^2 + 2 * (corNozzle*spotNozzle + divNozzle * z) * divNozzle * z - divNozzle^2*z^2);
    
    % assign energy to stf and run MC simulation
    stf.ray.energy = currentEnergy;
    
    %% needs to use correct BDL file in calcParticleDoseMC
    resultGUI = matRad_calcDoseDirectMC(ct,stf,pln,cst,1);          
    
    machine.data(count) = matRad_fitBaseData(resultGUI.physicalDose, ct.resolution, currentEnergy, mcData);
    
    matRad_plotParticleBaseDataEntry(machine,count,hf);
    count = count + 1;
end

%% save final machine  
matRad_cfg.dispInfo('Saving final machine %s to %s\n',fileName,filePath);
save(filePath,'machine','-v7');
clear machine; 
