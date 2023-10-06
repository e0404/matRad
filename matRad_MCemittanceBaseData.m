classdef matRad_MCemittanceBaseData
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % matRad_MCemmitanceBaseData This is the superclass for MonteCarlo base
    % data calculation
    %
    %
    %
    %
    %
    % References
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2020 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        machine                             %matRad base data machine struct
        bdl_path = ''                       %stores path to generated file
        nozzleToIso                         %Nozzle to Isocenter Distance
        smx                                 %Scanning magnet X to isocenter Distance
        smy                                 %Scanning magnet y to isocenter Distance
        monteCarloData                      %MC Phase space data struct
        selectedFocus                       %array containing selected focus indices per energy
        defaultRelativeEnergySpread = 0;    %default energy spread
        matRad_cfg                          %matRad config
        rangeShifters                       %Stores range shifters
        
        % air correction in beam optics approximation
        fitWithSpotSizeAirCorrection  = true;
        fitCorrectDoubleGaussian      = true;

        %To force the phase space approximation even if we have the data
        forceSpectrumApproximation  = false; 
        forceEmittanceApproximation = false;
        forceFixedMU                = true;
    end
    
    properties (SetAccess = private)
        stfCompressed   %measure whether function has additional info about
        %the stf
        problemSigma    % = 1, when there was a problem calculating sigma
        energyIndex     %Indices of calculated energies
    end
    
    methods
        function obj = matRad_MCemittanceBaseData(machine,stf)
            %matRad_MCsquareBaseData construct an instance of the MCsquare
            %Base data format using a focus index
            
            %stfCompressed states whether monteCarloData are calculated for
            %all energies (false) or only for energies which exist in given
            %stf. If function is called without stf stfCompressed = false.
            if nargin < 2 || isempty(stf)
                obj.stfCompressed = false;
            else
                obj.stfCompressed = true;
                obj = obj.getRangeShiftersFromStf(stf);
            end
            
            obj.matRad_cfg = MatRad_Config.instance();
                        
            obj.machine = machine;
            obj.problemSigma = false;
            obj.selectedFocus = ones(numel(machine.data),1) * NaN;
            
            if isfield(machine.meta,'BAMStoIsoDist')
                obj.nozzleToIso = machine.meta.BAMStoIsoDist;
            else
                obj.matRad_cfg.dispWarning('No information on BAMS to isocenter distance. Using generic value of 500mm');
                obj.nozzleToIso = 500;
                obj.machine.meta.BAMStoIsoDist = 500;
            end
            
            if all(isfield(machine.meta,{'SAD_x','SAD_y'}))
                obj.smx = machine.meta.SAD_x;
                obj.smy = machine.meta.SAD_y;
            elseif isfield(machine.meta,'SAD')
                SAD = machine.meta.SAD;
                obj.smx = SAD;
                obj.smy = SAD;
            else
                obj.matRad_cfg.dispError('No SAD found!');
            end
            
            obj.monteCarloData = [];
            
            %select needed energies and according focus indices by using stf
            if obj.stfCompressed
                tmp = [stf(:).ray];
                plannedEnergies     = [tmp.energy];
                focusIndex          = [tmp.focusIx];
                [~, ind]            = unique(plannedEnergies);
                plannedEnergies     = plannedEnergies(ind);
                focusIndex          = focusIndex(ind);
                [~ ,obj.energyIndex, ~] = intersect([machine.data(:).energy],plannedEnergies);
                
                %if no stf was refered all energies are chosen, while setting
                %the focus index for all energies to preliminary 1
            else
                plannedEnergies = [machine.data(:).energy];
                focusIndex = ones(size(plannedEnergies));
                [~ ,obj.energyIndex, ~] = intersect([machine.data(:).energy],plannedEnergies);
            end
            
            obj.selectedFocus(obj.energyIndex) = focusIndex;
            
            count = 1;
            for i = 1:length(obj.energyIndex)
                ixE = obj.energyIndex(i);
                
                %look up whether MonteCarlo data are already present in
                %machine file , if so do not recalculate
                if isfield(machine.data(ixE), 'energySpectrum') && ~obj.forceSpectrumApproximation
                    energySpectrum = machine.data(ixE).energySpectrum;
                    if isfield(energySpectrum,'type') && strcmp(energySpectrum.type,'gaussian')
                        energyData.NominalEnergy    = ones(1,4) * machine.data(ixE).energy(:);
                        energyData.MeanEnergy       = machine.data(ixE).energySpectrum.mean(:);
                        energyData.EnergySpread     = machine.data(ixE).energySpectrum.sigma(:);
                    else
                        energyData = obj.fitPhaseSpaceForEnergy(ixE);
                    end
                else
                    energyData = obj.fitPhaseSpaceForEnergy(ixE);
                end

                if isfield(machine.data(ixE),'MUdata') && ~obj.forceFixedMU
                    energyData.ProtonsMU = machine.data(ixE).MUdata.numParticlesPerMU;
                else
                    energyData.ProtonsMU = 1e6; %Standard w calibration
                end
                
                if isfield(machine.data(ixE).initFocus,'emittance') && ~obj.forceEmittanceApproximation 
                    data = [];
                    focusIx = obj.selectedFocus(ixE);
                    emittance = machine.data(ixE).initFocus.emittance(focusIx);

                    if ~strcmpi(emittance.type,'bigaussian')
                        matRad_cfg.dispError('Can not handle emittance of type ''%S''!',emittance.type);
                    end

                    nGauss = 1;
                    if isfield(emittance,'weight') 
                        nGauss = length(emittance.weight) + 1;
                    end

                    if nGauss > 2
                        matRad_cfg.dispError('Can not process more than two Gaussians in Emittance parameterization!');
                    end
                    
                    opticsData.Weight1          = 1;
                    opticsData.SpotSize1x       = emittance.sigmaX(1);
                    opticsData.Divergence1x     = emittance.divX(1);
                    opticsData.Correlation1x    = emittance.corrX(1);
                    opticsData.SpotSize1y       = emittance.sigmaY(1);
                    opticsData.Divergence1y     = emittance.divY(1);
                    opticsData.Correlation1y    = emittance.corrY(1);
                    
                    if nGauss == 1
                        opticsData.Weight2          = 0;
                        opticsData.SpotSize2x       = 0;
                        opticsData.Divergence2x     = 0;
                        opticsData.Correlation2x    = 0;
                        opticsData.SpotSize2y       = 0;
                        opticsData.Divergence2y     = 0;
                        opticsData.Correlation2y    = 0;
                    else
                        opticsData.Weight1          = 1 - emittance.weight(1);
                        opticsData.Weight2          = emittance.weight(1);
                        opticsData.SpotSize2x       = emittance.sigmaX(2);
                        opticsData.Divergence2x     = emittance.divX(2);
                        opticsData.Correlation2x    = emittance.corrX(2);
                        opticsData.SpotSize2y       = emittance.sigmaY(2);
                        opticsData.Divergence2y     = emittance.divY(2);
                        opticsData.Correlation2y    = emittance.corrY(2);
                    end

                    %opticsData.FWHMatIso = 2.355 * sigmaNull;
                    opticsData.FWHMatIso = machine.data(ixE).initFocus.SisFWHMAtIso;
                                
                    tmp = energyData;
                    f = fieldnames(opticsData);
                    for a = 1:length(f)
                        tmp.(f{a}) = opticsData.(f{a});
                    end
                    
                    data = [data; tmp];
                else
                    data = [];
                    tmp = energyData;
                    for j = 1:size(machine.data(ixE).initFocus.sigma,1)
                        
                        %                         tmp = energyData;
                        opticsData = obj.fitBeamOpticsForEnergy(ixE, j);
                        
                        f = fieldnames(opticsData);
                        for a = 1:length(f)
                            if j == 1
                                tmp.(f{a}) = opticsData.(f{a});
                            else
                                tmp.(f{a}) = [tmp.(f{a}), opticsData.(f{a})];
                            end
                        end
                        
                        data = tmp;
                    end
                    
                end

                                
                obj.monteCarloData = [obj.monteCarloData, data];
                count = count + 1;
            end
            
            %throw out warning if there was a problem in calculating the
            %width of the Bragg peak in obj.fitBeamOpticsForEnergy
            if obj.problemSigma
                obj.matRad_cfg.dispWarning('Calculation of FWHM of bragg peak in base data not possible! Using simple approximation for energy spread');
            end
        end
        
        
        function mcDataEnergy = fitPhaseSpaceForEnergy(obj,energyIx)
            %function to calculate mean energy and energy spread used by
            %mcSquare for given energy
            
            %Considers air distance from nozzle to phantom surface
            %used in the machine data. 0 means fitted to vacuum simulations
            %with surface at isocenter
            if ~isfield(obj.machine.meta, 'fitAirOffset')
                fitAirOffset = 0;
                %               warning('Could not find fitAirOffset. Using default value (no correction / fit in vacuum).');
            else
                fitAirOffset = obj.machine.meta.fitAirOffset;
            end
            dR = 0.0011 * (fitAirOffset);
            %dR = 0;
            
            i = energyIx;
            
            mcDataEnergy.NominalEnergy = ones(1, size(obj.machine.data(1).initFocus.dist,1)) * obj.machine.data(i).energy;
            
            newDepths = linspace(0,obj.machine.data(i).depths(end),numel(obj.machine.data(i).depths) * 100);
            newDose   = interp1(obj.machine.data(i).depths, obj.machine.data(i).Z, newDepths, 'spline');
            
            %find FWHM w50 of bragg peak and range of 80% does fall off
            [maxV, maxI] = max(newDose);
            [~, r80ind] = min(abs(newDose(maxI:end) - 0.8 * maxV));
            r80ind = r80ind - 1;
            r80 = interp1(newDose(maxI + r80ind - 1:maxI + r80ind + 1), ...
                newDepths(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);% ...
            % + obj.machine.data(i).offset + dR;
            
            %Correct r80 with dR
            r80 = r80 + dR + obj.machine.data(i).offset;
            
            [~, d50rInd] = min(abs(newDose(maxI:end) - 0.5 * maxV));
            d50rInd = d50rInd - 1;
            d50_r = interp1(newDose(maxI + d50rInd - 1:maxI + d50rInd + 1), ...
                newDepths(maxI + d50rInd - 1:maxI + d50rInd + 1), 0.5 * maxV);
            
            if (newDose(1) < 0.5 * maxV)
                [~, d50lInd] = min(abs(newDose(1:maxI) - 0.5*maxV));
                d50_l = interp1(newDose(d50lInd - 1:d50lInd + 1), ...
                    newDepths(d50lInd - 1:d50lInd + 1), 0.5 * maxV);
                w50 = d50_r - d50_l;
                %if width left of peak cannot be determined use r80 as width
            else
                % d50_l = newDepths(maxI);
                w50 = r80;
                obj.problemSigma = true;
            end
            
            %material dependent straggling factor for water (Bortfeld 1998
            %Eq. (18) & (B2)
            %We assume alpha to contain the A/Z^2 dependence across ions!
            alphaStraggling = 0.086; %MeV^2/cm
            fStragglingFactor = @(alpha,p) sqrt(alphaStraggling * p^3 * alpha^(2/p) / (3*p - 2)); %(sigma_straggling = 0.012*R^0.935), for calculations in cm
            fEnergyFromRange = @(R,alpha,p) (R./(10*alpha)).^(1/p);
            
                        
            %calcualte mean energy used my mcSquare with a formula fitted
            %to TOPAS data
            switch obj.machine.meta.radiationMode
                case 'protons'
                    %some constants according to Bortfeld (1993). Note that
                    %the paper usually works in [cm]
                    alpha = 2.2e-3;
                    p = 1.77;
                    
                    stragglingFactor = fStragglingFactor(alpha,p);
                    
                    %some functions describing range/energy relation and
                    %straggling

                    %polyfit to MC for energy from range (in mm)
                    fMeanEnergyFromRange = @(R) 5.762374661332111e-20 * R.^9 - 9.645413625310569e-17 * R.^8 + 7.073049219034644e-14 * R.^7 ...
                        - 2.992344292008054e-11 * R.^6 + 8.104111934547256e-09 * R.^5 - 1.477860913846939e-06 * R.^4 ...
                        + 1.873625800704108e-04 * R.^3 - 1.739424343114980e-02 * R.^2 + 1.743224692623838e+00 * R ...
                        + 1.827112816899668e+01;                 
                    %alternatively we can use range energy relationship
                    %Bortfeld (1993) Eq. (4) * 10 for [mm]
                    %meanEnergyFromRange = @(R) (R./(10*alpha)).^(1/p);

                    %polyfit to MC for straggling width
                    %sigmaRangeStragglingOnlySq = @(x) 2.713311945114106e-20 * x^9 - 4.267890251195303e-17 * x^8 + 2.879118523083018e-14 * x^7 ...
                    %     - 1.084418008735459e-11 * x^6 + 2.491796224784373e-09 * x^5 - 3.591462823163767e-07 * x^4 ...
                    %     + 3.232810400304542e-05 * x^3 - 1.584729282376364e-03 * x^2 + 5.228413840446568e-02 * x ...
                    %     - 6.547482267336220e-01;
                    %
                    
                    %straggling contribution according to Bortfeld Eq.
                    %(18), in [mm]
                    fStragglingSigmaFromRange = @(R) 10 * stragglingFactor * (R/10)^((3-2/p)/2);

                    %energy spectrum contribution to peak width according
                    %to Bortfeld Eq. 19, in mm
                    fEnergySpreadFromWidth = @(sigmaSq,E) sqrt(sigmaSq ./ ((10*alpha)^2 * p^2 * E^(2*p-2)));

                    mcDataEnergy.MeanEnergy = fMeanEnergyFromRange(r80);
                    
                    %calculate energy straggling using formulae deducted from paper
                    %"An analytical approximation of the Bragg curve for therapeutic
                    %proton beams" by T. Bortfeld et al. After inversion of
                    %the formula to obtain the two values z_50 where
                    %d(z_50) = 0.5*dMax, we obtain that the width is 6.14 *
                    %the total (energy + range) straggling sigma
                    totalSigmaSq = ((w50) / 6.14)^2;
                    
                    %Obtain estimate straggling component
                    sigmaRangeStragglingOnlySq = fStragglingSigmaFromRange(r80).^2; 

                    %Squared difference to obtain residual width from
                    %energy spectrum
                    if totalSigmaSq > sigmaRangeStragglingOnlySq
                        sigmaEnergyContributionSq = totalSigmaSq - sigmaRangeStragglingOnlySq;
                        energySpreadInMeV = fEnergySpreadFromWidth(sigmaEnergyContributionSq,mcDataEnergy.MeanEnergy);                
                    else
                        energySpreadInMeV = 1e-8; %monoenergetic, but let's not write 0 to avoid division by zero in some codes
                    end

                    energySpreadRelative = energySpreadInMeV ./ mcDataEnergy.MeanEnergy * 100;
                        
                    mcDataEnergy.EnergySpread = energySpreadRelative;
                case 'carbon'
                    %Constants
                    alpha = 4.425e-3;
                    p = 1.64;

                    alphaStraggling = 0.086; %MeV^2/cm
                    
                    % Fit to Range-Energy relationship
                    % Data from "Update to ESTAR, PSTAR, and ASTAR Databases" - ICRU Report 90, 2014
                    % Normalized energy before fit (MeV/u)! Only used ranges [10 350] mm for fit
                    % https://www.nist.gov/system/files/documents/2017/04/26/newstar.pdf
                    %meanEnergyFromRange = @(R) 11.39 * R^0.628 + 11.24;
                    fMeanEnergyFromRange = @(R) fEnergyFromRange(R,alpha,p);
                    mcDataEnergy.MeanEnergy = fMeanEnergyFromRange(r80);
                    % reading in a potential given energyspread could go here directly. How would you parse the energyspread
                    % into the function? Through a field in the machine?
                    
                    %Straggling factor:
                    %stragglingSigmaFromRange = @(R) 10 * stragglingFactor * (R/10)^0.935;


                    mcDataEnergy.EnergySpread = obj.defaultRelativeEnergySpread;
                case 'helium'
                    alpha = 2.567e-3;
                    p = 1.74;                                

                    % Fit to Range-Energy relationship
                    % Data from "Update to ESTAR, PSTAR, and ASTAR Databases" - ICRU Report 90, 2014
                    % Normalized energy before fit (MeV/u)! Only used ranges [10 350] mm for fit
                    % https://www.nist.gov/system/files/documents/2017/04/26/newstar.pdf
                    %meanEnergyFromRange = @(x) 7.57* x.^0.5848 + 3.063;

                    fMeanEnergyFromRange = @(R) fEnergyFromRange(R,alpha,p);
                    mcDataEnergy.MeanEnergy = fMeanEnergyFromRange(r80);


                    stragglingFactor = fStragglingFactor(alpha,p);
                    fStragglingSigmaFromRange = @(R) 10 * stragglingFactor * (R/10)^((3-2/p)/2);
                    fEnergySpreadFromWidth = @(sigmaSq,E) sqrt(sigmaSq ./ ((10*alpha)^2 * p^2 * E^(2*p-2)));
                    totalSigmaSq = ((w50) / 6.14)^2;
                    sigmaRangeStragglingOnlySq = fStragglingSigmaFromRange(r80).^2; 

                    if totalSigmaSq > sigmaRangeStragglingOnlySq
                        sigmaEnergyContributionSq = totalSigmaSq - sigmaRangeStragglingOnlySq;
                        energySpreadInMeV = fEnergySpreadFromWidth(sigmaEnergyContributionSq,mcDataEnergy.MeanEnergy);                
                    else
                        energySpreadInMeV = 1e-8; %monoenergetic, but let's not write 0 to avoid division by zero in some codes
                    end

                    energySpreadRelative = energySpreadInMeV ./ mcDataEnergy.MeanEnergy * 100;
                        
                    mcDataEnergy.EnergySpread = energySpreadRelative;

                    %mcDataEnergy.EnergySpread = obj.defaultRelativeEnergySpread;
                otherwise
                    error('not implemented')
            end
        end
        
        
        
        
        function mcDataOptics = fitBeamOpticsForEnergy(obj,energyIx, focusIndex)
            %function to calculate beam optics used by mcSquare for given
            %energy
            
            i = energyIx;
            
            %calculate geometric distances and extrapolate spot size at nozzle
            SAD = obj.machine.meta.SAD;
            z     = -(obj.machine.data(i).initFocus.dist(focusIndex,:) - SAD);
            sigma = obj.machine.data(i).initFocus.sigma(focusIndex,:);
            
            %Double Gaussian data might have a non-zero wide Gaussian,
            %adding width to the beam. We do a maximum correction here,
            %which compromises the fwhm, but seems to work better in
            %estimating the optics
            if obj.fitCorrectDoubleGaussian && isfield(obj.machine.data(i),'sigma1')              
                sigmaSq_Narr = sigma.^2 + obj.machine.data(i).sigma1(1).^2;
                sigmaSq_Bro  = sigma.^2 + obj.machine.data(i).sigma2(1).^2;
                dgWeight = obj.machine.data(i).weight(1);
                
                %Maximum of double gaussian
                maxL = (1-dgWeight) ./ (2*pi*sigmaSq_Narr) + dgWeight ./ (2*pi*sigmaSq_Bro );
                
                %Find the sigma that corresponds to the maximum
                sigma = sqrt(1 ./ (2*pi*maxL));
            end
            
            %correct for in-air scattering with polynomial or interpolation
            if obj.fitWithSpotSizeAirCorrection
                sigma = arrayfun(@(d,sigma) obj.spotSizeAirCorrection(obj.machine.meta.radiationMode,obj.machine.data(i).energy,d,sigma),-z+obj.machine.meta.BAMStoIsoDist,sigma);                   
            end

            %square and interpolate at isocenter
            sigmaSq = sigma.^2;     
            sigmaSqIso = sqrt(interp1(z,sigmaSq,0));
            
            %fit Courant-Synder equation to data using ipopt, formulae
            %given in mcSquare documentation            
            
            %fit function
            qRes = @(rho, sigmaT) (sigmaSq -  (sigmaSqIso^2 - 2*sigmaSqIso*rho*sigmaT.*z + sigmaT^2.*z.^2));

            % fitting for either matlab or octave
            if ~obj.matRad_cfg.isOctave
                funcs.objective = @(x) sum(qRes(x(1), x(2)).^2);
                funcs.gradient  = @(x) [  2 * sum(qRes(x(1), x(2)) .* (2 * sigmaSqIso * x(2) * z));
                    2 * sum(qRes(x(1), x(2)) .* (2 * sigmaSqIso * x(1) * z  - 2 * x(2) * z.^2))];
                
                options.lb = [-0.99, -Inf];
                options.ub = [ 0.99,  Inf];
                
                options.ipopt.hessian_approximation = 'limited-memory';
                options.ipopt.limited_memory_update_type = 'bfgs';
                options.ipopt.print_level = 1;
                
                start = [0.9; 0.1];
                [result, ~] = ipopt (start, funcs, options);
                rho    = result(1);
                sigmaT = result(2);
                
            else
                phi{1} = @(x) sum(qRes(x(1), x(2)).^2);
                phi{2} = @(x) [  2 * sum(qRes(x(1), x(2)) .* (2 * sigmaSqIso * x(2) * z));
                    2 * sum(qRes(x(1), x(2)) .* (2 * sigmaSqIso * x(1) * z  - 2 * x(2) * z.^2))];
                
                lb = [-0.99, -Inf];
                ub = [ 0.99,  Inf];
                
                start = [0.9; 0.1];
                [result, ~] = sqp (start, phi, [], [], lb, ub);
                rho    = result(1);
                sigmaT = result(2);
            end
            
            %calculate divergence, spotsize and correlation at nozzle
            DivergenceAtNozzle  = sigmaT;
            SpotsizeAtNozzle    = sqrt(sigmaSqIso^2 - 2 * rho * sigmaSqIso * sigmaT * obj.nozzleToIso + sigmaT^2 * obj.nozzleToIso^2);
            CorrelationAtNozzle = (rho * sigmaSqIso - sigmaT * obj.nozzleToIso) / SpotsizeAtNozzle;
            
            
            %save calcuated beam optics data in mcData
            mcDataOptics.ProtonsMU     = 1e6;
            
            mcDataOptics.Weight1       = 1;
            mcDataOptics.SpotSize1x    = SpotsizeAtNozzle;
            mcDataOptics.Divergence1x  = DivergenceAtNozzle;
            mcDataOptics.Correlation1x = CorrelationAtNozzle;
            mcDataOptics.SpotSize1y    = SpotsizeAtNozzle;
            mcDataOptics.Divergence1y  = DivergenceAtNozzle;
            mcDataOptics.Correlation1y = CorrelationAtNozzle;
            
            visBool = false;
            if visBool
                figure, plot(z,sigmaSq,'x');
                zNew = linspace(z(1),z(end),100);
                y = sigmaSqIso^2 - 2*rho*sigmaSqIso*sigmaT * zNew + sigmaT^2 * zNew.^2;
                hold on; plot(zNew,y);
            end
            
            
            mcDataOptics.Weight2       = 0;
            mcDataOptics.SpotSize2x    = 0;
            mcDataOptics.Divergence2x  = 0;
            mcDataOptics.Correlation2x = 0;
            mcDataOptics.SpotSize2y    = 0;
            mcDataOptics.Divergence2y  = 0;
            mcDataOptics.Correlation2y = 0;
            mcDataOptics.FWHMatIso = 2.355 * sigmaSqIso;
        end
        
        
        
        function obj = saveMatradMachine(obj,name)
            %save previously calculated monteCarloData in new baseData file
            %with given name
            
            %             [~ ,energyIndex, ~] = intersect([obj.machine.data(:).energy], [obj.monteCarloData(:).NominalEnergy]);
            
            machineName = [obj.machine.meta.radiationMode, '_', name];
            
            count = 1;
            for i = 1:length(obj.energyIndex)
                
                ixE = obj.energyIndex(i);

                numFoci = numel(obj.monteCarloData(:,count).SpotSize1x);

                for f = 1:numFoci
                    if obj.monteCarloData(:,count).Weight2 ~= 0
                        obj.machine.data(ixE).initFocus.emittance(f).type  = 'bigaussian';
                        obj.machine.data(ixE).initFocus.emittance(f).sigmaX = [obj.monteCarloData(f,count).SpotSize1x obj.monteCarloData(:,count).SpotSize2x(f)];
                        obj.machine.data(ixE).initFocus.emittance(f).sigmaY = [obj.monteCarloData(f,count).SpotSize1y obj.monteCarloData(:,count).SpotSize2y(f)];
                        obj.machine.data(ixE).initFocus.emittance(f).divX = [obj.monteCarloData(f,count).Divergence1x obj.monteCarloData(:,count).Divergence2x(f)];
                        obj.machine.data(ixE).initFocus.emittance(f).divY = [obj.monteCarloData(f,count).Divergence1y obj.monteCarloData(:,count).Divergence2y(f)];
                        obj.machine.data(ixE).initFocus.emittance(f).corrX = [obj.monteCarloData(f,count).Correlation1x obj.monteCarloData(:,count).Correlation2x(f)];
                        obj.machine.data(ixE).initFocus.emittance(f).corrY = [obj.monteCarloData(f,count).Correlation1y obj.monteCarloData(:,count).Correlation2y(f)];
                        obj.machine.data(ixE).initFocus.emittance(f).weight = [obj.monteCarloData(f,count).Weight2]; %Weight one will not be stored explicitly due to normalization
                    else
                        obj.machine.data(ixE).initFocus.emittance(f).type  = 'bigaussian';
                        obj.machine.data(ixE).initFocus.emittance(f).sigmaX = [obj.monteCarloData(:,count).SpotSize1x(f)];
                        obj.machine.data(ixE).initFocus.emittance(f).sigmaY = [obj.monteCarloData(:,count).SpotSize1y(f)];
                        obj.machine.data(ixE).initFocus.emittance(f).divX = [obj.monteCarloData(:,count).Divergence1x(f)];
                        obj.machine.data(ixE).initFocus.emittance(f).divY = [obj.monteCarloData(:,count).Divergence1y(f)];
                        obj.machine.data(ixE).initFocus.emittance(f).corrX = [obj.monteCarloData(:,count).Correlation1x(f)];
                        obj.machine.data(ixE).initFocus.emittance(f).corrY = [obj.monteCarloData(:,count).Correlation1y(f)];
                    end
                end
                
                %At the moment coded to only take the first energy because
                %focus settings do not apply to the energy spectrum
                obj.machine.data(ixE).energySpectrum.type  = 'gaussian';
                obj.machine.data(ixE).energySpectrum.mean   = [obj.monteCarloData(:,count).MeanEnergy];
                obj.machine.data(ixE).energySpectrum.sigma = [obj.monteCarloData(:,count).EnergySpread];
                
                
                count = count + 1;
            end
            machine = obj.machine;
            machineFilePath = fullfile(obj.matRad_cfg.matRadRoot,'basedata',[machineName '.mat']);

            save('-v7',machineFilePath,'machine');
            obj.matRad_cfg.dispInfo('Saved Emittance to matRad base data in %s\n',machineFilePath);
        end
    end
    
    methods (Access = protected)
        function obj = getRangeShiftersFromStf(obj,stf)
            allRays = [stf.ray];
            raShis = [allRays.rangeShifter];
            
            [~,ix] =  unique(cell2mat(squeeze(struct2cell(raShis))'),'rows');
            
            raShis = raShis(ix);
            
            ix = [raShis.ID] == 0;
            
            obj.rangeShifters = raShis(~ix);
        end
    end

    methods (Static)
        function sigmaAirCorrected = spotSizeAirCorrection(radiationMode,E,d,sigma,method)
            %performs a rudimentary correction for additional scattering in
            %air not considered by the courant snyder equation
            matRad_cfg = MatRad_Config.instance();
            if nargin < 5
                method = 'interp_linear';
            end

            switch radiationMode
                case 'protons'
                    %Provide Look-up Table and fit for protons
                    sigmaLUT = [0    0.4581    2.7777    7.0684   12.6747; ...
                                0    0.1105    0.7232    2.1119    4.2218; ...
                                0    0.0754    0.5049    1.4151    2.8604; ...
                                0    0.0638    0.3926    1.1196    2.2981; ...
                                0    0.0466    0.3279    0.9440    1.9305; ...
                                0    0.0414    0.2825    0.8294    1.7142; ...
                                0    0.0381    0.2474    0.7336    1.5192; ...
                                0    0.0335    0.2214    0.6696    1.3795; ...
                                0    0.0287    0.2030    0.6018    1.2594; ...
                                0    0.0280    0.1925    0.5674    1.1865; ...
                                0    0.0257    0.1801    0.5314    1.0970; ...
                                0    0.0244    0.1670    0.4966    1.0342];
                    energies = [31.7290   69.4389   95.2605  116.5270  135.1460  151.9670  167.4620  181.9230  195.5480  208.4780  220.8170  232.6480]';
                    depths = [0 500 1000 1500 2000];                                       
                    polyFit = @(E,d) 0.001681*d - 0.0001178*E*d + 6.094e-6*d^2 + 1.764e-6*E^2*d - 1.016e-7*E*d^2 - 9.803e-09*E^3*d + 6.096e-10*E^2*d^2 + 1.835e-11*E^4*d - 1.209e-12*E^3*d^2;
                otherwise 
                    %No air correction because we don't have data yet
                    sigmaLUT = [0 0; 0 0];
                    energies = [0; realmax];
                    depths = [0; realmax];

                    polyFit = @(E,d) 0;
            end

            %make sure to not violate ranges!
            %this is a little hardcoded, but helps us handle strange
            %distances in the initFocus field
            if d < min(depths)
                d = min(depths);
                matRad_cfg.dispWarning('Spot Size Air Correction problem, negative distance found!',method);
            end

            if d > max(depths)
                d = max(depths);
                matRad_cfg.dispWarning('Spot Size Air Correction problem, distance too large!',method);
            end

            if E > max(energies)
                E = max(energies);
                matRad_cfg.dispWarning('Spot Size Air Correction problem, energy too large!',method);
            end

            if E < min(energies)
                E = min(energies);
                matRad_cfg.dispWarning('Spot Size Air Correction problem, energy too small!',method);
            end


            switch method
                case 'interp_linear'
                    sigmaAir = interp2(energies,depths,sigmaLUT',E,d,'linear');
                case 'fit'
                    sigmaAir = polyFit(E,d);
                otherwise
                    matRad_cfg.dispWarning('Air Correction Method ''%s'' not known, skipping!',method);
                    sigmaAir = 0;
            end
    
            if sigmaAir >= sigma
                sigmaAirCorrected = sigma;
                matRad_cfg.dispWarning('Spot Size Air Correction failed, too large!',method);
            else
                sigmaAirCorrected = sqrt(sigma.^2 - sigmaAir.^2); 
            end

        end
    end
end

