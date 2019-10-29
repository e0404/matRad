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
load BOXPHANTOMv2.mat

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
%pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propStf.isoCenter       = [ct.cubeDim(1) / 2 * ct.resolution.x, ...
                                0, ...%ct.cubeDim(2) / 2 * ct.resolution.y, ...
                                ct.cubeDim(3) / 2 * ct.resolution.z];

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = ct.resolution.x*3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = ct.resolution.y*3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = ct.resolution.z*3; % [mm]

% optimization settings
pln.propOpt.optimizer       = 'IPOPT';
pln.propOpt.bioOptimization = 'none'; % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                      % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;  % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

%% generate steering file
stf = matRad_generateStf(ct,cst,pln);
stf.ray.energy = 110;      %50 bis 220  

 %% dose calculation
Ecount = 1;
for i = 50:50:150
    stf.ray.energy = i
    dijMC = matRad_calcParticleDoseMC(ct,stf,pln,cst);
    resultGUI = matRad_calcCubes(ones(dijMC.totalNumOfBixels,1),dijMC);           

    IDD = sum(sum(resultGUI.physicalDose,2),3);
    IDD = reshape(IDD,320,1);
    IDD = [IDD(1); IDD];
    depthsIDD = ct.resolution.y/2 : ct.resolution.y : ct.resolution.y*ct.cubeDim(2);


    depthsGauss = depthsIDD - (depthsIDD(end)-depthsIDD(1))/2;

    depthsIDD = [0, depthsIDD];

    % plot(depthsIDD, IDD);

    % figure;
    % plot(depthsGauss,profile);

    %% curve fitting
    resSigma = [];
    for i = 1:320
        profile = resultGUI.physicalDose(i,:,160);

        gauss = @(x, a, sigma) a * exp(- x.^2 / (2*sigma^2));

        funcs.objective = @(p) sum((gauss(depthsGauss, p(1), p(2)) - profile).^2);

        funcs.gradient = @(p) [ sum(2 * (p(1) * exp(-depthsGauss.^2 / (2 * p(2)^2)) - profile) ...
                                                .* exp(-depthsGauss.^2 / (2 * p(2)^2)));
                                sum(2 * (p(1) * exp(-depthsGauss.^2 / (2 * p(2)^2)) - profile) ...
                                                .* p(1) .* exp(-depthsGauss.^2 / (2 * p(2)^2)) .* depthsGauss.^2 / p(2)^3)];

        options.lb = [0,  0];
        options.ub = [ Inf,  Inf];

        options.ipopt.hessian_approximation = 'limited-memory';
        options.ipopt.limited_memory_update_type = 'bfgs';
        options.ipopt.print_level = 1;

        start = [1; 1];
        [result, ~] = ipopt (start, funcs, options);


    %     plot(depthsGauss,gauss(depthsGauss,result(1),result(2)));
    %     hold on
    %     scatter(depthsGauss,profile);
    %     waitforbuttonpress
        resSigma = [resSigma, result(2)];
    end   



    resSigma = [resSigma(1), resSigma];
    initSigma = sum(resSigma(1:10)) / 10;
    resSigma = sqrt(resSigma.^2 - initSigma^2);
    resSigma(imag(resSigma) > 0) = 0;              



    %interpolate range at 80% dose after peak.
                [maxV, maxI] = max(IDD);
                [~, r80ind] = min(abs(IDD(maxI:end) - 0.8 * maxV));
                r80ind = r80ind - 1;
                r80 = interp1(IDD(maxI + r80ind - 1:maxI + r80ind + 1), ...
                                 depthsIDD(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);


    machine.meta.radiationMode = 'protons';
    machine.meta.dataType = 'singleGauss';
    machine.meta.created_on = date;
    machine.meta.created_by = 'pamede';
    machine.meta.SAD = (2218 + 1839) / 2;
    machine.meta.BAMStoIsoDist = 420.0;
    machine.meta.machine = 'Generic';
    machine.meta.LUT_bxWidthminFWHM = [1, Inf; initSigma, initSigma];


    machine.data(Ecount).range = r80;    

    machine.data(Ecount).energy =  stf.ray.energy;

    machine.data(Ecount).depths = depthsIDD';

    machine.data(Ecount).Z = IDD;

    machine.data(Ecount).peakPos = depthsIDD(maxI);

    machine.data(Ecount).sigma = resSigma';

    machine.data(Ecount).offset = 0;

    machine.data(Ecount).initFocus.dist  = [0, 20000];
    machine.data(Ecount).initFocus.sigma = [initSigma, initSigma];
    machine.data(Ecount).initFocus.SisFWHMAtIso = 2.3548 * initSigma;

    Ecount = Ecount + 1;
end
% nozzleToIso = 420.0;
% 
% SAD = (2218 + 1839) / 2;
% 
% z = nozzleToIso;
% 
% SpotSize1x    = 5.154355331;  
% Divergence1x  = 0.004403330;
% Correlation1x = 0.457081065;
% 
% 
% spotsizeIso = sqrt(SpotSize1x^2 + 2 * Divergence1x * z * ...
%     (Correlation1x * SpotSize1x - Divergence1x * z) - Divergence1x^2 * z^2);
%                 
% resSigma = sqrt(resSigma.^2 - spotsizeIso^2);
% 
% plot(resSigma)













% matRadGUI;
