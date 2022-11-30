classdef MatRad_HeterogeneityConfig < handle
    % MatRad_TopasConfig class definition
    %
    %
    % References
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2019 the matRad development team.
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
        % Instance of matRad configuration class
        matRad_cfg = MatRad_Config.instance();

        % Heterogeneity correction can be turned off with setting this flag to 'false'
        calcHetero = true;

        useOriginalDepths = false;
        bioOpt = false;

        modulateLET = true;
        modulateBioDose = true;

        % "medium" modulation power
        % Pmod = 256; % [µm]
        % worst case modulation power
        modPower = 800;

        useDoseCurves = true;

        type = 'complete';  % 'complete','depthBased','voxelwise'
        
        % Property struct for sampling
        sampling = struct('mode','matRad', ...
            'method','binomial',...  % 'binomial','poisson'
            ...%         'numOfHistories',1e6,...
            'numOfSamples',20,...
            'continuous',true);

    end

    methods
        function obj = MatRad_HeterogeneityConfig()
            % MatRad_TopasConfig Construct configuration Class for TOPAS
        end

        function out = Gauss(~,x,mu,SqSigma)
            % Function handle for calculating lateral dose
            out = 1./(sqrt(2*pi.*SqSigma)).*exp(-((x - mu).^2./(2.*SqSigma)));
        end

        function out = sumGauss(~,x,mu,SqSigma,w)
            % Function handle for calculating depth doses
            out = (1./sqrt(2*pi*ones(numel(x),1) .* SqSigma') .* exp(-bsxfun(@minus,x,mu').^2 ./ (2* ones(numel(x),1) .* SqSigma' ))) * w;
        end

        function sigmaSq = getHeterogeneityCorrSigmaSq(obj,WET,Pmod)
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % matRad calculation of Bragg peak degradation due to heterogeneities
            %
            % call
            %   sigmaSq = obj.getHeterogeneityCorrSigmaSq(WET)
            %
            % input
            %   WET:        water equivalent thickness of heterogeneous structure, e.g. lung [mm]
            %   Pmod:       modulation power [µm] (optional)
            %
            % output
            %   sigmaSq:    sigma squared [mm^2] for degradation of Bragg peak
            %
            % References
            %   -
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % a = 1.43 [mm/sqrt(cm)] from Riccardos thesis, fit MC well,
            % independent from proton energy, R80 stays at same position with degradation
            % a = 1.60 [mm/sqrt(cm)] = 1.6/sqrt(10) [sqrt(mm)] fits the measurements better
            %
            % Pmod = a^2; 150-750 micrometer for swine lung from Witt et al.
            % Pmod = 256 [µm] is equivalent to a = 1.6 [mm/sqrt(cm)]
            % Pmod = 204.5 [µm] to a = 1.43 [mm/sqrt(cm)]
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Copyright 2019 the matRad development team.
            %
            % This file is part of the matRad project. It is subject to the license
            % terms in the LICENSE file found in the top-level directory of this
            % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
            % of the matRad project, including this file, may be copied, modified,
            % propagated, or distributed except according to the terms contained in the
            % LICENSE file.
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if nargin < 3
                Pmod = obj.modPower;
            end
            % output sigma^2 in mm^2
            sigmaSq = Pmod/1000 .* WET;
        end

        function cst = cstHeteroAutoassign(~,cst)
            % Prepares the cst file for the heterogeneity correction algorithm
            %
            % call
            %   cstHetero = obj.cstHeteroAutoassign(cst)
            %
            % input
            %   cst:      matRad cst struct
            %
            % output
            %   cst:      updated matRad cst struct with 'Lung' property
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Copyright 2019 the matRad development team.
            %
            % This file is part of the matRad project. It is subject to the license
            % terms in the LICENSE file found in the top-level directory of this
            % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
            % of the matRad project, including this file, may be copied, modified,
            % propagated, or distributed except according to the terms contained in the
            % LICENSE file.
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Note: tumor tissue does not contribute to the degradation, but some seams around the GTV might.
            % lungTissue={'Lung','GTV','PTV','CTV','ITV'};
            lungTissue={'Lung'};

            % assign the 'Lung' property to the segmentations containing the string "lung".
            for i = 1:length(cst(:,1))
                if any(cellfun(@(teststr) ~isempty(strfind(cst{i,2},teststr)), lungTissue))
                    cst{i,5}.HeterogeneityCorrection = 'Lung';
                end
            end

        end

        function samples = sampleBino(obj,n,p,numOfSamples,continuous)
            % matRad function for binomial sampling
            %
            % call
            %   X = obj.sampleBino(n,p,numOfSamples)
            %
            % input
            %   n:              number of independent experiments
            %   p:              probability (between 0 and 1)
            %   numOfSamples:   number of samples for output
            %
            % output
            %   X:              binomial samples
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

            % handle inputs
            if any(~mod(n,1) == 0) && ~continuous
                obj.matRad_cfg.dispError('n has to be integer for discrete distribution');
            elseif ~(all(p <= 1) && all(p >= 0))
                obj.matRad_cfg.dispError('p must be between 0 and 1');
            elseif ~isscalar(p) && ~(numOfSamples == numel(p))
                obj.matRad_cfg.dispError('p array must have numOfSamples entries')
            elseif ~isscalar(n) && ~(numOfSamples == numel(n))
                obj.matRad_cfg.dispError('n array must have numOfSamples entries')
            end

            % save time when testing on homogeneous phantom or when performing density
            % override
            if isscalar(unique(n)) && isscalar(unique(p))
                n = unique(n);
                p = unique(p);
            end

            if continuous
                % calculate beta distribution parameters alpha & beta using method of moments with sample mean x=p and
                % https://en.wikipedia.org/wiki/Beta_distribution#Method_of_moments
                a = p .* (n-1);
                b = (1-p) .* (n-1);

                %     a = p.*n;
                %     b = (1-p).*n;

                % sample from continuous beta distribution using "numOfSamples" random numbers
                samples = betaincinv(rand([numOfSamples,1]),a,b);
            else
                % sample discrete binomial distribution
                if isscalar(n)
                    samples = sum(rand([numOfSamples,n]) < p, 2);
                else
                    samples = zeros(numOfSamples,1);
                    for i = 1:numel(n)
                        samples(i) = sum(rand([1,n(i)]) < p(i), 2);
                    end
                end
                % need normalization here to only get values over [0,1]
                samples = samples ./ n;
            end

        end

        function ct = modulateDensity(obj,ct,cst,pln)
            % matRad density modulation function to calculate sampled ct struct
            %
            % call
            %   ct = obj.modulateDensity(ct,cst,Pmod,mode)
            %
            % input
            %   ct:             ct struct
            %   cst:            matRad cst struct
            %   Pmod:           Modulation power according to which the modulation will
            %                   be created
            %   mode:           mode for density modulation ('binominal','poisson')
            %                   note: poisson only available for Pmod = 250,450 and 800
            %
            % output
            %   ct:             ct struct with modulated density cube
            %
            % References
            %   [1] Poisson sampling: https://iopscience.iop.org/article/10.1088/1361-6560/aa641f
            %   [2] Detwiler et al., 2021, Compendium of Material Composition Data for Radiation Transport Modeling
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

            % get all unique lung indices from lung segmentations
            idx = cellfun(@(teststr) ~isempty(strfind(lower(teststr),'lung')), cst(:,2));
            if sum(idx)==0
                obj.matRad_cfg.dispError('No lung segmentation found in cst.\n');
            end
            lungIdx = [cst{idx,4}];
            lungIdx = unique(vertcat(lungIdx{:}));

            % calculate ct cube from cubeHU if not specified
            if ~isfield(ct,'cube')
                ct = matRad_calcWaterEqD(ct,pln);
            end

            switch pln.propHeterogeneity.sampling.method
                case 'binomial'
                    % set lung tissue density
                    rhoLung = 1.05; % [2] [g/cm^3]

                    % isolate lung densities and normalize to lung tissue density
                    pLung = ct.cube{1}(lungIdx) / rhoLung;

                    % scrap lung densities larger than 1 (cannot be used in sampling, since distributions are defined over [0,1])
                    if any(pLung > 1)
                        lungIdx = lungIdx(pLung <= 1);
                        pLung = ct.cube{1}(lungIdx) / rhoLung;
                    end

                    % calculate size of substructures
                    d = pln.propHeterogeneity.modPower/1000 ./ (1-pLung) / rhoLung; % [1] eq.8: Pmod = d*(1-pLung) * rhoLung

                    % length of a voxel (uniform voxels are assumed)
                    D = ct.resolution.y;

                    % calculate number of substructures inside of a single voxel (round in case of discrete distribution)
                    if pln.propHeterogeneity.sampling.continuous
                        n = D./d;
                        largeEnough = n > 1;
                    else
                        n = round(D./d);
                        largeEnough = n >= 1;
                    end

                    % Don't modulate voxel with less than 1 substructures
                    lungIdx = lungIdx(largeEnough);
                    pLung = pLung(largeEnough);
                    n = n(largeEnough);

                    % get samples from the binomial distribution (discrete or continuous approximation)
                    samples = obj.sampleBino(n,pLung,length(lungIdx),pln.propHeterogeneity.sampling.continuous);

                    % revert normalization to get values between [0,rhoLung]
                    samples = samples * rhoLung;

                    % write samples to CT and convert to Hounsfield Units
                    ct.cube{1}(lungIdx) = samples;
                    ct.cubeHU{1}(lungIdx) = 1024*(ct.cube{1}(lungIdx)-1);

                case 'poisson'
                    % read density modulation look-up table for supported modulation powers (Pmod = 250, 450, 800 mu)
                    DensMod = obj.loadModDist(pln.propHeterogeneity.modPower);

                    % Connect each density-probability-couple with a number that will later be transformed the HU-Value:
                    % The maximum HU-Value of the HU set by Schneider etal is 2995: So the HU of the modulated density must be at least 2995+1=2996 ; This values is prerocessed with the later used RescaleIntercept and RescaleSlope. See also calculation for Threshold_HU_Value_to_double.
                    Min_HU_for_DensMod_in_double=((2995+1000)/1);

                    % Real HU values in row 3
                    DensMod(:,3) = (Min_HU_for_DensMod_in_double+1:Min_HU_for_DensMod_in_double+length(DensMod))' - 1000;

                    % Use 50 times as many samples since samples are being scrapped in the next step
                    numOfSamples = 50 * length(lungIdx);

                    % Sample density modulation
                    zz1 = ceil(rand(numOfSamples,1)*length(DensMod));
                    zz2 = rand(numOfSamples,1);
                    ix = zz2 <= DensMod(zz1,2);
                    newDist = zz1(ix);

                    % Write samples in HU to CT
                    ct.cubeHU{1}(lungIdx) = DensMod(newDist(1:numel(lungIdx)),3);

                    % Write density cube in CT
                    %     ct.cube{1}(lungIdx) = (ct.cubeHU{1}(lungIdx)+1)/1024;

                    % Descrete sampling of the density distribution
                    %     P = [0; cumsum(DensMod(:,2))];
                    %     samples = discretize(rand(numel(lungIdx),1),P);
                    %     ct.cube{1}(lungIdx) = samples / max(samples);
            end

            % In case of Monte Carlo simulations, the Hounsfield Units have to be adjusted to work with the respective
            % density to material converters.
            if any(cellfun(@(teststr) ~isempty(strfind(pln.propHeterogeneity.sampling.mode,teststr)), {'TOPAS','MCsquare'}))
                % Only include different densities that are significantly different (1e-3)
                % This is done to significantly increase the computation time
                lung = round(ct.cube{1}(lungIdx),3);

                % Set minimum air density to 0.001225 to avoid deviding by 0
                lung(lung < 0.001225) = 0.001225;

                % Extract distinct densities and respective number of occurences
                [numOfOccurences,sampledDensities] = groupcounts(lung);

                %
                switch pln.propMC.materialConverter.addSection
                    % Modes used to include samples in material converter in case of discrete sampling (either 0 or 1.05)
                    case 'lung'
                        ct.cubeHU{1}(lungIdx(lung == 1.05)) = 0;
                        ct.cubeHU{1}(lungIdx(lung == 0))    = -999;
                    case 'poisson'
                        ct.cubeHU{1}(lungIdx(lung == 1.05)) = 3020;
                        ct.cubeHU{1}(lungIdx(lung == 0))    = 2997;
                        % Main mode used by TOPAS and MCsquare pipeline to include samples in material converter
                    case 'sampledDensities'
                        % Sort densities
                        [~,sortIdx] = sort(lung);

                        % Set individual "virtual" Hounsfield Units (beginning at 6000) for each density and repeat for their number of occurence
                        lungDensitiesNewSorted = repelem(6000:5999+numel(numOfOccurences),1,numOfOccurences);

                        % Write sorted lung densities to HU cube (using also sorted indices)
                        ct.cubeHU{1}(lungIdx(sortIdx)) = lungDensitiesNewSorted;

                        % Save sampled densities and indices to CT to write in material converter
                        ct.sampledDensities   = sampledDensities;
                        ct.sampledLungIndices = lungIdx(sortIdx);
                    otherwise
                        obj.matRad_cfg.dispWarning('Lung modulation should be used with a separate section in the Schneider converter.\n');
                end
            end

            % Set flag to indicate that the CT has been modulated
            ct.modulated = 1;

            % plot histogram of the the lung density distribution
            % figure, histogram(ct.cube{1}(lungIdx))
        end

        function resultGUI = accumulateOverSamples(~,resultGUI,resultGUI_mod,samples)
            % Don't accumulate RBE by default
            accumulateRBE = false;

            % Get number of beams from resultGUI_mod
            fnames = fieldnames(resultGUI_mod);
            
            fnames = fnames(cellfun(@(teststr) ~isempty(strfind(teststr,'beam')), fnames));
            fnames = cellfun(@(x) strsplit(x,'_'), fnames, 'UniformOutput', false);
            fnames = [fnames{:,1}];
            fnames = erase(unique(fnames(cellfun(@(teststr) ~isempty(strfind(teststr,'beam')), fnames))),'beam');
            numOfBeams = max(cellfun(@(x) str2double(x), fnames));

            % get beam info
            for i = 1:numOfBeams
                beamInfo(i).suffix = ['_beam', num2str(i)];
            end
            beamInfo(numOfBeams+1).suffix = '';

            % Load RBE models if MonteCarlo was calculated for multiple models
            if any(cellfun(@(teststr) ~isempty(strfind(teststr,'RBExD')), fieldnames(resultGUI_mod)))
                accumulateRBE = true;
            end

            % Handle RBE models if available
            if accumulateRBE
                if isfield(resultGUI_mod,'RBE_model')
                    RBE_model = cell(1,length(resultGUI_mod.RBE_model)+1);
                    for i = 1:length(resultGUI_mod.RBE_model)
                        RBE_model{i+1} = ['_' resultGUI_mod.RBE_model{i}];
                    end
                else
                    RBE_model = {''};
                end
            end

            % Allocate empty dose fields if resultGUI is empty
            if isempty(fieldnames(resultGUI))
                for i = 1:length(beamInfo)
                    resultGUI.(['physicalDose' beamInfo(i).suffix]) = zeros(size(resultGUI_mod.physicalDose));

                    if accumulateRBE
                        for j = 1:length(RBE_model)
                            resultGUI.(['RBExD' RBE_model{j} beamInfo(i).suffix]) = zeros(size(resultGUI_mod.physicalDose));
                        end
                    end
                end
            end

            % Accumulate averaged physical Dose
            for i = 1:length(beamInfo)
                resultGUI.(['physicalDose' beamInfo(i).suffix]) = resultGUI.(['physicalDose' beamInfo(i).suffix]) + resultGUI_mod.(['physicalDose' beamInfo(i).suffix]) / samples;
            end

            % Accumulate averaged RBE weighted Dose
            if accumulateRBE
                for i = 1:length(beamInfo)
                    for j = 1:length(RBE_model)
                        resultGUI.(['RBExD' RBE_model{j} beamInfo(i).suffix]) = resultGUI.(['RBExD' RBE_model{j} beamInfo(i).suffix]) + resultGUI_mod.(['RBExD' RBE_model{j} beamInfo(i).suffix]) / samples;
                    end
                end
            end

        end

        function stdOut = calcSampleStd(~,dataPoints,mean)

            meanDiff = 0;
            samples = length(dataPoints);

            for k = 1:samples
                meanDiff = meanDiff + (dataPoints{k} - mean).^2;
            end
            varMean = meanDiff./(samples - 1)./samples;
            stdMean = sqrt(varMean);

            stdSum = stdMean * samples;
            varSum = stdSum.^2;

            stdOut = sqrt(varSum);

        end
    end

    methods (Access = private)

        function DensMod = loadModDist(~,Pmod)
            % matRad function to load density modulation tables
            % (only for Pmod = 250, 450 and 800)
            %
            % call
            %   [DensMod] = matRad_loadModDist(Pmod)
            %
            % input
            %   Pmod:           Modulation Power
            %
            % output
            %   DensMod:        Density probabilty distribution
            %
            % References
            %   [1] https://iopscience.iop.org/article/10.1088/1361-6560/aa641f
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

            DensMod = load(['DensMod_',num2str(Pmod),'mu_0,2603_rho_1,5mm_pixelsize.txt']);
            DensMod(:,2) = DensMod(:,2) / sum(DensMod(:,2));

        end
    end

    methods (Static)

        function baseData = overrideBaseData(baseData)
            % script to switch from APM (machine.data.Z is a struct) to "original"
            % calculation (machine.data.Z is double)
            %
            % call
            %   baseData = obj.overrideBaseData(baseData)
            %
            % input
            %   baseData:       matRad machine base data
            %
            % output
            %   baseData        updated base data with new depths if possible
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Copyright 2019 the matRad development team.
            %
            % This file is part of the matRad project. It is subject to the license
            % terms in the LICENSE file found in the top-level directory of this
            % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
            % of the matRad project, including this file, may be copied, modified,
            % propagated, or distributed except according to the terms contained in the
            % LICENSE file.
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            matRad_cfg = MatRad_Config.instance();

            if isstruct(baseData(1).Z) && isfield(baseData(1).Z,'profileORG')
                matRad_cfg.dispInfo('Overriding fitted APM data with original depth dose.\n');
                for i = 1:length(baseData)
                    baseData(i).Z = baseData(i).Z.profileORG;
                end
            elseif isstruct(baseData(1).Z) && isfield(baseData(1).Z,'doseORG')
                matRad_cfg.dispInfo('Overriding fitted APM data with original depth dose.\n');
                for i = 1:length(baseData)
                    baseData(i).Z = baseData(i).Z.doseORG;
                end
            elseif isstruct(baseData(1).Z) && (~isfield(baseData(1).Z,'profileORG') && ~isfield(baseData(1).Z,'doseORG'))
                warning('No original depths available in base data. Nothing changed.');
            else
                warning('Base data depths are already in the desired format.');
            end

        end
    end
end

