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
% load PHANTOM_all.mat

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

stf = matRad_generateStf(ct,cst,pln);
 
 %% Monte Carlo dose
load protons_matRadBDL_APM;
stf.ray.energy = machine.data(42).energy;                           

resultGUI_MC = matRad_calcDoseDirectMC(ct,stf,pln,cst,ones(sum(stf(:).totalNumOfBixels),1), 5000000);
resultGUI.physicalDoseMC = resultGUI_MC.physicalDose;
mcDose      = resultGUI.physicalDoseMC;  
  
figure
imagesc(mcDose(:,:,round(ct.cubeDim(3)/2)));
caxis([0 2e-3]);
title('Monte Carlo dose');
hold on 
contour(ct.cube{1}(:,:,round(ct.cubeDim(3)/2)),2,'color','white');
contour(mcDose(:,:,round(ct.cubeDim(3)/2)),10,'color','black');
hold off
pbaspect([ct.cubeDim(2) ct.cubeDim(1) ct.cubeDim(3)])

 %% 
 
gamma = [];
sigma = [];
for i = 4.5:0.1:6
    pln.propDoseCalc.anaMode = 'stdCorr';
    dijSC = matRad_calcParticleDose(ct,stf,pln,cst,false, i);
    resultGUI_SC = matRad_calcCubes(ones(sum(stf(:).totalNumOfBixels),1),dijSC);
    resultGUI.physicalDoseSC = resultGUI_SC.physicalDose;
    anaScDose   = resultGUI.physicalDoseSC;

    gammaTest = [2, 2];
    interpolation = 0;
    [~,gammaPassRateCell] = matRad_gammaIndex(anaScDose,mcDose,[ct.resolution.x,ct.resolution.y,ct.resolution.z],gammaTest,round(ct.cubeDim(3)/2),interpolation,'global',cst);
    gamma = [gamma, gammaPassRateCell{1,2}];
    sigma = [sigma, i]
end
    
plot(sigma,gamma)  

