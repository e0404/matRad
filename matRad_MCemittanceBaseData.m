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
        rangeShifters                       %Stores range shifters

        % air correction in beam optics approximation
        fitWithSpotSizeAirCorrection  = true;
        fitCorrectDoubleGaussian      = false;

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

            % Instance of MatRad Config
            matRad_cfg = MatRad_Config.instance();

            %stfCompressed states whether monteCarloData are calculated for
            %all energies (false) or only for energies which exist in given
            %stf. If function is called without stf stfCompressed = false.
            if nargin < 2 || isempty(stf)
                obj.stfCompressed = false;
            else
                obj.stfCompressed = true;
                obj = obj.getRangeShiftersFromStf(stf);
            end

            obj.machine = machine;
            obj.problemSigma = false;
            obj.selectedFocus = ones(numel(machine.data),1) * NaN;

            if isfield(machine.meta,'BAMStoIsoDist')
                obj.nozzleToIso = machine.meta.BAMStoIsoDist;
            else
                matRad_cfg.dispWarning('No information on BAMS to isocenter distance. Using generic value of 500mm');
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
                matRad_cfg.dispError('No SAD found!');
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

            % Loop through all required energies
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

                %% Skip emittance approximation if emittance was found in base data
                if isfield(machine.data(ixE).initFocus,'emittance') && ~obj.forceEmittanceApproximation
                    data = [];
                    focusIx = obj.selectedFocus(ixE);
                    emittance = machine.data(ixE).initFocus.emittance(focusIx);

                    if ~strcmpi(emittance.type,'bigaussian')
                        matRad_cfg.dispError('Can not handle emittance of type ''%S''!',emittance.type);
                    end

                    if isfield(emittance,'weight')
                        nGauss = length(emittance.weight) + 1;
                    else
                        nGauss = 1;
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
            end

            %throw out warning if there was a problem in calculating the
            %width of the Bragg peak in obj.fitBeamOpticsForEnergy
            if obj.problemSigma
                matRad_cfg.dispWarning('Calculation of FWHM of bragg peak in base data not possible! Using simple approximation for energy spread');
            end
        end


        function mcDataEnergy = fitPhaseSpaceForEnergy(obj,energyIx)
            % function to calculate mean energy and energy spread used by Monte Carlo for given energy

            % Instance of MatRad Config
            matRad_cfg = MatRad_Config.instance();

            % Considers air distance from nozzle to phantom surface used in the machine data.
            % 0 means fitted to vacuum simulations with surface at isocenter
            if ~isfield(obj.machine.meta, 'fitAirOffset')
                fitAirOffset = 0;
                % warning('Could not find fitAirOffset. Using default value (no correction / fit in vacuum).');
            else
                fitAirOffset = obj.machine.meta.fitAirOffset;
            end
            % Introduce air offset correction
            airOffsetCorrection = 0.0011 * (fitAirOffset);

            % Save given energy as nominal energy
            mcDataEnergy.NominalEnergy = obj.machine.data(energyIx).energy;

            % Find range of 80% does fall off after the peak
            [maxDose, maxDoseIdx] = max(obj.machine.data(energyIx).Z);
            % interpolation to evaluate interpolated depths at 80% maxDose (constrain interpolation to area after peak)
            r80 = interp1(obj.machine.data(energyIx).Z(maxDoseIdx:end), obj.machine.data(energyIx).depths(maxDoseIdx:end), 0.8 * maxDose);
            % Correct r80 with air offset and potential offset from basedata
            r80 = r80 + airOffsetCorrection + obj.machine.data(energyIx).offset;

            % Calcualte mean energy used my mcSquare with a formula fitted to TOPAS data
            switch obj.machine.meta.radiationMode
                case 'protons'
                    %%% Approximate mean energy
                    % This rangeEnergy relationship was created analogously to the helium and carbon relationship below
                    % Bragg-Kleeman rule R(E_0)=\alpha E_0^p (Bragg 1905) inversely fitted to data.
                    % Only used ranges [10 300] mm for fit. Units: [E]=MeV/u, [R]=mm.
                    % data from Berger2023 https://dx.doi.org/10.18434/T4NC7P
                    meanEnergyFromRange = @(R) 8.531* R.^0.5663;

                    %%% Calculate energy spread from FWHM
                    % Calculate FWHM of bragg peak if the plateau is lower than 50% of the max Dose
                    if (obj.machine.data(energyIx).Z(1) < 0.5 * maxDose)
                        % Calculate left and right flank, where dose falls to 50%
                        d50_r = interp1(obj.machine.data(energyIx).Z(maxDoseIdx:end), obj.machine.data(energyIx).depths(maxDoseIdx:end), 0.5 * maxDose);
                        d50_l = interp1(obj.machine.data(energyIx).Z(1:maxDoseIdx), obj.machine.data(energyIx).depths(1:maxDoseIdx), 0.5 * maxDose);

                        % Combine for FWHM
                        FWHM = d50_r - d50_l;
                    else
                        % if width left of peak cannot be determined use r80 as width
                        FWHM = r80;
                        obj.problemSigma = true;
                        matRad_cfg.warning('Could not find FWHM, using r80 as approximation');
                    end

                    % Calculate energy straggling using formulae deducted from paper
                    % "An analytical approximation of the Bragg curve for therapeutic proton beams" by T. Bortfeld et al.
                    % After inversion of the formula to obtain the two values z_50 where d(z_50) = 0.5*dMax,
                    % we obtain that the width is 6.14 * the total (energy + range) straggling sigma
                    totalSigmaSq = (FWHM / 6.14)^2;

                    % Bortfeld 1997, Eq.17
                    % Changed alpha and p based on new energy-range fit made analogously to range-energy fit above
                    % (original values are commented behind)
                    alpha   = 0.02383;  % alpha = 0.022; [alpha] = mm MeV^-p
                    p       = 1.756;    % p = 1.77; [p] = 1
                    alphaPrime = 0.0087; % (MeV^2/mm)

                    sigmaRangeStragglingOnlySq = @(R) alphaPrime * (p^3*alpha^(2/p))/(3*p-2) * R ^(3-2/p);

                    % Use formula deducted from Bragg Kleeman rule to calcuate energy straggling given the total sigma
                    % and the range straggling (Bortfeld 1997, Eq.19, in mm)
                    energySpreadFromWidth = @(sigmaSq,E) sqrt(sigmaSq ./ (alpha^2 * p^2 * E^(2*p-2)));

                    %Squared difference to obtain residual width from energy spectrum
                    if totalSigmaSq > sigmaRangeStragglingOnlySq(r80)
                        sigmaEnergyContributionSq = totalSigmaSq - sigmaRangeStragglingOnlySq(r80);
                        energySpreadInMeV = energySpreadFromWidth(sigmaEnergyContributionSq,mcDataEnergy.MeanEnergy);
                    else
                        energySpreadInMeV = 1e-8; %monoenergetic, but let's not write 0 to avoid division by zero in some codes
                    end

                    energySpreadRelative = energySpreadInMeV ./ mcDataEnergy.MeanEnergy * 100;
                case 'carbon'
                    %%% Approximate mean energy
                    % Fit to Range-Energy relationship
                    % Data from "Update to ESTAR, PSTAR, and ASTAR Databases" - ICRU Report 90, 2014
                    % Normalized energy before fit (MeV/u)! Only used ranges [10 300] mm for fit
                    % https://www.nist.gov/system/files/documents/2017/04/26/newstar.pdf
                    meanEnergyFromRange = @(R) 14.35 * R^0.5915;

                    %                     alpha = 0.0127;
                    %                     p = 1.667;

                    %%% Energy spread
                    % Reading in a potential given energyspread could go here directly. How would you parse the energyspread
                    % into the function? Through a field in the machine?
                    % TODO this could be done analogously to protons, adjusting the alpha and p parameter from a E0(R0)
                    % graph (inverse of the fit above)
                    energySpreadRelative = obj.defaultRelativeEnergySpread;

                case 'helium'
                    %%% Fit to Range-Energy relationship
                    % Data from "Update to ESTAR, PSTAR, and ASTAR Databases" - ICRU Report 90, 2014
                    % Normalized energy before fit (MeV/u)! Only used ranges [10 300] mm for fit
                    % https://www.nist.gov/system/files/documents/2017/04/26/newstar.pdf
                    meanEnergyFromRange = @(R) 8.465* R.^0.5672;

                    %                     alpha = 0.02461;
                    %                     p = 1.751;

                    %%% Energy spread
                    energySpreadRelative = obj.defaultRelativeEnergySpread;

                otherwise
                    error('not implemented')
            end

            % Write previously approximated meanEnergy and energySpread
            mcDataEnergy.MeanEnergy = meanEnergyFromRange(r80);
            mcDataEnergy.EnergySpread = energySpreadRelative;
        end




        function mcDataOptics = fitBeamOpticsForEnergy(obj,energyIx, focusIndex)
            % function to calculate beam optics used by Monte Carlo for given energy

            % Instance of MatRad Config
            matRad_cfg = MatRad_Config.instance();

            %calculate geometric distances and extrapolate spot size at nozzle
            SAD = obj.machine.meta.SAD;
            z     = -(obj.machine.data(energyIx).initFocus.dist(focusIndex,:) - SAD);
            sigma = obj.machine.data(energyIx).initFocus.sigma(focusIndex,:);

            % TODO Add comment
            if obj.fitCorrectDoubleGaussian
                sigmaSq_Narr = sigma.^2 + obj.machine.data(i).sigma1(1).^2;
                sigmaSq_Bro  = sigma.^2 + obj.machine.data(i).sigma2(1).^2;
                sigma = sqrt(sigmaSq_Narr*(1-obj.machine.data(i).weight(1).^2) + sigmaSq_Bro*obj.machine.data(i).weight(1).^2);
            end

            %correct for in-air scattering with polynomial or interpolation
            sigma = arrayfun(@(d,sigma) obj.spotSizeAirCorrection(obj.machine.meta.radiationMode,obj.machine.data(energyIx).energy,d,sigma),-z+obj.machine.meta.BAMStoIsoDist,sigma);

            %square and interpolate at isocenter
            sigmaSq = sigma.^2;
            sigmaSqIso = sqrt(interp1(z,sigmaSq,0));

            %fit Courant-Synder equation to data using ipopt, formulae
            %given in mcSquare documentation

            %fit function
            qRes = @(rho, sigmaT) (sigmaSq -  (sigmaSqIso^2 - 2*sigmaSqIso*rho*sigmaT.*z + sigmaT^2.*z.^2));

            % Define optimization parameters

            start = [0.9; 0.1];
            options.lb = [-0.99, -Inf];
            options.ub = [ 0.99,  Inf];

            funcs.objective = @(x) sum(qRes(x(1), x(2)).^2);
            funcs.gradient  = @(x) [  2 * sum(qRes(x(1), x(2)) .* (2 * sigmaSqIso * x(2) * z));
                2 * sum(qRes(x(1), x(2)) .* (2 * sigmaSqIso * x(1) * z  - 2 * x(2) * z.^2))];

            % fitting for either matlab or octave
            if ~matRad_cfg.isOctave
                options.ipopt.hessian_approximation = 'limited-memory';
                options.ipopt.limited_memory_update_type = 'bfgs';
                options.ipopt.print_level = 1;

                [result, ~] = ipopt (start, funcs, options);
            else
                [result, ~] = sqp (start, phi, [], [], lb, ub);

            end

            rho    = result(1);
            sigmaT = result(2);

            %calculate divergence, spotsize and correlation at nozzle
            DivergenceAtNozzle  = sigmaT;
            SpotsizeAtNozzle    = sqrt(sigmaSqIso^2 - 2 * rho * sigmaSqIso * sigmaT * obj.nozzleToIso + sigmaT^2 * obj.nozzleToIso^2);
            CorrelationAtNozzle = (rho * sigmaSqIso - sigmaT * obj.nozzleToIso) / SpotsizeAtNozzle;


            % Save calcuated beam optics data in mcData
            mcDataOptics.ProtonsMU     = 1e6;

            mcDataOptics.Weight1       = 1;
            mcDataOptics.SpotSize1x    = SpotsizeAtNozzle;
            mcDataOptics.Divergence1x  = DivergenceAtNozzle;
            mcDataOptics.Correlation1x = CorrelationAtNozzle;
            mcDataOptics.SpotSize1y    = SpotsizeAtNozzle;
            mcDataOptics.Divergence1y  = DivergenceAtNozzle;
            mcDataOptics.Correlation1y = CorrelationAtNozzle;
            mcDataOptics.Weight2       = 0;
            mcDataOptics.SpotSize2x    = 0;
            mcDataOptics.Divergence2x  = 0;
            mcDataOptics.Correlation2x = 0;
            mcDataOptics.SpotSize2y    = 0;
            mcDataOptics.Divergence2y  = 0;
            mcDataOptics.Correlation2y = 0;

            visBool = false;
            if visBool
                figure, plot(z,sigmaSq,'x');
                zNew = linspace(z(1),z(end),100);
                y = sigmaInit^2 - 2*rho*sigmaInit*sigmaT * zNew + sigmaT^2 * zNew.^2;
                hold on; plot(zNew,y);
            end

            mcDataOptics.FWHMatIso = 2.355 * sigmaSqIso;
        end



        function obj = saveMatradMachine(obj,name)
            %save previously calculated monteCarloData in new baseData file
            %with given name

            % Instance of MatRad Config
            matRad_cfg = MatRad_Config.instance();

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
                obj.machine.data(ixE).energySpectrum.mean   = [obj.monteCarloData(:,count).MeanEnergy(f)];
                obj.machine.data(ixE).energySpectrum.sigma = [obj.monteCarloData(:,count).EnergySpread(f)];


                count = count + 1;
            end
            machine = obj.machine;
            machineFilePath = fullfile(matRad_cfg.matRadRoot,'basedata',[machineName '.mat']);

            save('-v7',machineFilePath,'machine');
            matRad_cfg.dispInfo('Saved Emittance to matRad base data in %s\n',machineFilePath);
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

            % Instance of MatRad Config
            matRad_cfg = MatRad_Config.instance();

            if nargin < 5
                method = 'interp_linear';
            end

            switch radiationMode
                case 'protons'
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
                    sigmaLUT = [0 0; 0 0];
                    energies = [0; Inf];
                    depths = [0; Inf];

                    polyFit = @(E,d) 0;
            end

            %make sure to not violate ranges!
            %this is a little hardcoded, but helps us handle strange distances in the initFocus field
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
                sigmaAirCorrected = sigma - sigmaAir;
            end

        end
    end
end

