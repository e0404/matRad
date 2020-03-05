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
% load PHANTOM_slab_entrance_10mm.mat
load PHANTOM_slab_entrance_50mm.mat
% load PHANTOM_slab_entrance_100mm.mat
% load PHANTOM_fineSamplingTest.mat
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
energyIx = 37;
stf.ray.energy = machine.data(energyIx).energy;
                            
%% dose calculation
 % analytical dose without fine sampling
    tic
    pln.propDoseCalc.anaMode = 'standard';
    dij = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI = matRad_calcCubes(ones(sum(stf(:).totalNumOfBixels),1),dij);
    anaDose     = resultGUI.physicalDose;
    t1 = toc;

 % Monte Carlo dose
    tic
    resultGUI_MC = matRad_calcDoseDirectMC(ct,stf,pln,cst,ones(sum(stf(:).totalNumOfBixels),1), 1000000);
    resultGUI.physicalDoseMC = resultGUI_MC.physicalDose;
    mcDose      = resultGUI.physicalDoseMC;
    t4 = toc;

%% radDepths fitting
%%    
meanRadDepths(meanRadDepths > machine.data(energyIx).depths(end)) = machine.data(energyIx).depths(end);
ergRadDepth = zeros(ct.cubeDim);
ergSigma = zeros(ct.cubeDim);
for i = 80:140%1:ct.cubeDim(1)
    i
    for j = 1:ct.cubeDim(2)
        for k = 25%1:ct.cubeDim(3)
            f1 = @(x) abs(mcDose(i,j,k) - matRad_calcParticleDoseBixel(x(1), radialDist(i,j,k)^2, sigmaIni_sq, machine.data(energyIx), x(2))) * 100000000;
%             f2 = @(x) (mcDose(i,j,k) - matRad_calcParticleDoseBixel(meanRadDepths(i,j,k), radialDist(i,j,k)^2, sigmaIni_sq, machine.data(energyIx), x(1)))^2 * 100000000;
%             f3 = @(x) (mcDose(i,j,k) - matRad_calcParticleDoseBixel(x(1), radialDist(i,j,k)^2, sigmaIni_sq, machine.data(energyIx), 0))^2 * 100000000;
            
            if meanRadDepths(i,j,k) > machine.data(energyIx).depths(end)
                guessX = [machine.data(energyIx).depths(end), 3]
            else
                guessX = [meanRadDepths(i,j,k), 3];
            end
%             range = 10;
%             if meanRadDepths(i,j,k) < range
%                 lb = [0, 0];
%             else
%                 lb = [meanRadDepths(i,j,k) - range, 0];
%             end
%             
%             if meanRadDepths(i,j,k) > (machine.data(energyIx).depths(end) - range)
%                 ub = [machine.data(energyIx).depths(end), 10];
%             else
%                 ub = [meanRadDepths(i,j,k) + range, 10];
%             end
            
            lb = [50, 0];
            ub = [machine.data(energyIx).depths(end), 5];
               
            options = optimset('Display', 'none');
            tmp = fmincon(f1, guessX, [], [], [], [], lb, ub, [], options);
%             tmp = simulannealbnd(f1, guessX, lb, ub, options);
        
            ergRadDepth(i,j,k) = tmp(1);
            ergSigma(i,j,k) = tmp(2);

        end
    end
end

sigmaMapSlice = ergSigma(:,:,25);
% meanRadDepthSlice = meanRadDepths(:,:,25);
radialDistSlice = radialDist(:,:,25);
dose = matRad_calcParticleDoseBixel(reshape(ergRadDepth(:,:,25),150*50,1), reshape(radialDistSlice,150*50,1).^2, sigmaIni_sq, machine.data(energyIx), reshape(ergSigma(:,:,25),150*50,1));%reshape(sigma(:,:,25),150*50,1));
dose = reshape(dose,150,50);

sigmaMapSlice = imgaussfilt(sigmaMapSlice,1);
figure
imagesc(sigmaMapSlice)
hold on
contour(mcDose(:,:,round(ct.cubeDim(3)/2)),linspace(0,2e-3,10),'color','black');

figure
imagesc(ergRadDepth(:,:,25))
hold on
contour(mcDose(:,:,round(ct.cubeDim(3)/2)),linspace(0,2e-3,10),'color','black');

figure
imagesc(dose)
hold on
contour(mcDose(:,:,round(ct.cubeDim(3)/2)),linspace(0,2e-3,10),'color','black');












