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
  
minEnergy = 220;
maxEnergy = 225;
nEnergy   = 3;  

count = 1;
for currentEnergy = linspace(minEnergy, maxEnergy, nEnergy)
    %assign energy to stf and run MC simulation
    stf.ray.energy = currentEnergy;
    resultGUI = matRad_calcDoseDirectMC(ct,stf,pln,cst,ones(sum(stf(:).totalNumOfBixels),1),1000000);          

    %extract IDD and calculate according depths/lengths for IDD and Gaussian fit
    IDD = sum(sum(resultGUI.physicalDose,2),3);
    IDD = reshape(IDD,1100,1);
    IDD = [IDD(1); IDD];
    IDDnotZero = find(IDD);
    IDD = IDD(IDDnotZero);
    
    depthsIDD = 0 : ct.resolution.x : ct.resolution.x * ct.cubeDim(1);
    depthsIDD = depthsIDD(IDDnotZero);
    IDD = interp1(depthsIDD, IDD, 0:0.05:depthsIDD(end), 'spline');
    depthsIDD = 0:0.05:depthsIDD(end);
    
    axisGauss = linspace( -ct.resolution.y * ct.cubeDim(2) / 2, ct.resolution.y * ct.cubeDim(2) / 2, ct.cubeDim(2)) + ct.resolution.y/2;
    
    %calculate sigma/spotsize at isocenter using MC phase space data
    divNozzle = interp1(energyMC,divMC,currentEnergy); 
    corNozzle = interp1(energyMC,corMC,currentEnergy); 
    spotNozzle = interp1(energyMC,spotMC,currentEnergy); 
    z = 420;
    spotIso = sqrt(spotNozzle^2 + 2 * (corNozzle*spotNozzle + divNozzle * z) * divNozzle * z - divNozzle^2*z^2);
    
    resSigma = [];
    for i = 1:(numel(IDDnotZero) - 1)
    %curve fitting gaussians
        
        
        %save cast perpendicular to beam axis in profile
        profile = resultGUI.physicalDose(i,:,125);

        gauss = @(x, a, sigma) a * exp(- x.^2 / (2*sigma^2));

        
        %define fit parameters
        funcs.objective = @(p) sum((gauss(axisGauss, p(1), p(2)) - profile).^2);
        
        funcs.gradient = @(p) [ sum(2 * (p(1) * exp(-axisGauss.^2 / (2 * p(2)^2)) - profile) ...
                                                .* exp(-axisGauss.^2 / (2 * p(2)^2)));
                                sum(2 * (p(1) * exp(-axisGauss.^2 / (2 * p(2)^2)) - profile) ...
                                                .* p(1) .* exp(-axisGauss.^2 / (2 * p(2)^2)) .* axisGauss.^2 / p(2)^3)];

        options.lb = [0,  0];
        options.ub = [ Inf,  Inf];
    
        options.ipopt.hessian_approximation = 'limited-memory';
        options.ipopt.limited_memory_update_type = 'bfgs';
        options.ipopt.print_level = 1;
        options.ipopt.tol = 1e-16;
        
        %run fit and calculate actual sigma by squared substracting initial
        %sigma / spotsize
        ixMid = round(numel(profile)/2);
        start = [sum(profile(ixMid-1:ixMid+1))/3; spotIso];
        [fitResult, ~] = ipopt (start, funcs, options);
        
%         plot(axisGauss, gauss(axisGauss, fitResult(1), fitResult(2)))
%         hold on
%         plot(axisGauss, profile)
%         hold off

        if(fitResult(2) > spotIso)
            actualSigma = sqrt(fitResult(2)^2 - spotIso^2);
        else
            actualSigma = 0;
        end
        
        resSigma = [resSigma, actualSigma];
    end   

    depthsSigma = 0 : ct.resolution.x : ct.resolution.x * ct.cubeDim(1);
    depthsSigma = depthsSigma(IDDnotZero);
    
    resSigma = [0, resSigma];
    resSigma = interp1(depthsSigma, resSigma, depthsIDD);
    
    %interpolate range at 80% dose after peak.
    [maxV, maxI] = max(IDD);
    [~, r80ind] = min(abs(IDD(maxI:end) - 0.8 * maxV));
    r80ind = r80ind - 1;
    r80 = interp1(IDD(maxI + r80ind - 1:maxI + r80ind + 1), ...
                     depthsIDD(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);

    %conversion factor for IDDs
    cf = 1 / 1.6021766208e-02;
    
    %save data in machine
    machine.meta.radiationMode = 'protons';
    machine.meta.dataType = 'singleGauss';
    machine.meta.created_on = date;
    machine.meta.created_by = 'Paul Anton Meder';
    machine.meta.SAD = (2218 + 1839) / 2;
    machine.meta.BAMStoIsoDist = 420.0;
    machine.meta.machine = 'Generic';
    machine.meta.LUT_bxWidthminFWHM = [1, Inf; 2.3548 * spotIso, 2.3548 * spotIso];


    machine.data(count).range = r80;    

    machine.data(count).energy =  stf.ray.energy;

    machine.data(count).depths = depthsIDD';

    machine.data(count).Z = IDD * cf;

    machine.data(count).peakPos = depthsIDD(maxI);

    machine.data(count).sigma = resSigma';

    machine.data(count).offset = 0;

    machine.data(count).initFocus.dist  = [0, 20000];
    machine.data(count).initFocus.sigma = [spotIso, spotIso];
    machine.data(count).initFocus.SisFWHMAtIso = 2.3548 * spotIso;
    
    display(['baseData Progress :', ' ', num2str(round(count/nEnergy*100)), '%']);
    count = count + 1;
    toc
end

%save new machine
save('protons_testMachineFit', 'machine');