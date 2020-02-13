%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Expample script for baseData fitting to mcSquare simulation %%%%%%%

 %% 1) CT & CST creation
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
ct.cube{1} = ones(cubeDim) * 1;
ct.cube{1}(1,1,1) = 0; 

ct.resolution.x = 0.32;
ct.resolution.y = 0.32;
ct.resolution.z = 0.32;

ct.cubeDim = cubeDim;

ct.numOfCtScen  = 1;


%% Create a cubic phantom
iso = [800,125,125];

% create an ct image series with zeros
ct.cubeHU{1} = ones(ct.cubeDim) * 0;
ct.cubeHU{1}(1,1,1) = -1000; 

ct.hlut = [1,0;0,-1024];

% create body of full phantom size
mask = ones(ct.cubeDim);
cst{1,4}{1} = find(mask == 1);

%create target
centerP_corr = iso;
height_corr = 10;
width_corr = 10;
depth_corr = 10;


mask = zeros(ct.cubeDim);

for i=-height_corr/2+1:height_corr/2
    for j=-width_corr/2:width_corr/2
        for k=-depth_corr/2:depth_corr/2
            mask(centerP_corr(1)+i,centerP_corr(2)+j,centerP_corr(3)+k) = 1;
        end
    end
end
cst{2,4}{1} = find(mask == 1);

disp('CT creation done!');
clearvars -except ct cst


 %% 2) MCsquare computation and baseData fitting

% meta information for treatment plan
pln.radiationMode   = 'protons';     % either photons / protons / carbon
pln.machine         = 'matRadBDLold';

pln.numOfFractions  = 30;
 
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

%% generate steering file
stf = matRad_generateStf(ct,cst,pln);
 
%read MC phase space data
dataMC = readtable('BDL_matrad.txt');             
dataMC = dataMC{11:end,:};
tmp = [];
for i = 1:size(dataMC,1)    
    tmp = [tmp; strsplit(dataMC{i})];
end
energyMC    = str2double(tmp(:,1));
spotMC      = (str2double(tmp(:,6)) + str2double(tmp(:,9)))  / 2;
divMC       = (str2double(tmp(:,7)) + str2double(tmp(:,10))) / 2;
corMC       = (str2double(tmp(:,8)) + str2double(tmp(:,11))) / 2;
  
% define energy range
minEnergy = 70;
maxEnergy = 225;
nEnergy   = 88;  

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
    
    %assign energy to stf and run MC simulation
    stf.ray.energy = currentEnergy;
    
    %% needs to use correct BDL file in calcParticleDoseMC
    resultGUI = matRad_calcDoseDirectMC(ct,stf,pln,cst,ones(sum(stf(:).totalNumOfBixels),1),1000000);          
    
    machine.data(count) = matRad_fitBaseData(resultGUI.physicalDose, ct.resolution, currentEnergy, mcData);

    disp(['baseData Progress :', ' ', num2str(round(count/nEnergy*100)), '%']);
    count = count + 1;
end

% save data in machine
machine.meta.radiationMode = 'protons';
machine.meta.dataType = 'singleGauss';
machine.meta.created_on = date;
machine.meta.created_by = 'Paul Anton Meder';
machine.meta.SAD = (2218 + 1839) / 2;
machine.meta.BAMStoIsoDist = 420.0;
machine.meta.machine = 'Generic';
machine.meta.LUT_bxWidthminFWHM = [1, Inf; 8 ,8];
machine.meta.fitAirOffset = 420.0;
  
% save machine
save('protons_fitMachine.mat', 'machine');