classdef matRad_PhotonCollapsedConeEngine < DoseEngines.matRad_PencilBeamEngineAbstract
    % matRad_PhotonCollapsedConeEngine: photon dose calculation with
    %   collapsed cone convolution/superposition of a polyenergetic energy
    %   deposition kernel (Ahnesjo 1989)
    %
    %   For each bixel the TERMA is computed on the dose grid voxels returned
    %   by the ray tracer and transported with a set of collapsed cones whose
    %   polar axis is tilted along the beam axis or the bixel's ray. The transport is performed
    %   by the mex file matRad_collapsedConeMex or, if it is not compiled, by
    %   the MATLAB implementation matRad_collapsedConeTransport. The absolute
    %   dose is normalized to the analytical SVD pencil beam depth dose of the
    %   machine at a reference depth in water.
    %
    %   Configurable properties (pln.propDoseCalc):
    %
    %   - numPolarAngles: number of polar cone bins (divisor of 48)
    %   - numAzimuthalAngles: number of azimuthal cones per polar bin
    %   - scatterCutOff: lateral margin [mm] of the transport sub-volume around the TERMA of each bixel; dose is stored
    %     within max(geometricLateralCutOff, scatterCutOff) around each ray
    %   - normalizationDepth: depth [mm] of the absolute dose normalization
    %   - kernelFile: kernel data file (default: 6 MV water kernel)
    %   - useMex: use the compiled transport if available
    %   - numThreads: number of OpenMP threads of the mex file (0 = default)
    %   - hardeningCorrection: use the depth dependent effective attenuation mu0 / (1 + beta d) of the TERMA, with mu0
    %     and beta fitted so that the broad field depth dose matches the machine data
    %   - numScatterExponentials: number of exponential terms of the scatter kernel (1 or 2; the stored kernel is
    %     refitted if it differs)
    %   - kernelTilt: polar axis of the cones, 'beam' (beam axis, same transport lattice for all bixels of a beam) or
    %     'ray' (ray of each bixel, which causes lattice aliasing in the summed dose)
    %   - enableDijSampling: enable the lateral dij sampling
    %   - dijSampling: sampling parameters, see matRad_sampleDij
    %
    %   Limitations: the spectral hardening of the TERMA is approximated by a
    %   depth dependent effective attenuation coefficient (no off-axis
    %   softening, no spectral kernel) and the kernel is computed in water
    %   and scaled with the relative electron density.
    %
    % References
    %   [1] Ahnesjo, Med. Phys. 16, 577 (1989)
    %   [2] Bortfeld, Schlegel, Rhein, Med. Phys. 20, 311 (1993)
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2026 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Constant)
        possibleRadiationModes = {'photons'} % radiation modes the engine can process
        name = 'Collapsed Cone Convolution' % user readable name of the engine
        shortName = 'CCC' % identifier used in pln.propDoseCalc.engine
    end

    properties (SetAccess = public, GetAccess = public)
        numPolarAngles = 12 % number of polar cone bins, must divide the 48 cones of the kernel data
        numAzimuthalAngles = 8 % number of azimuthal cones per polar bin
        scatterCutOff = 150 % lateral margin [mm] of the transport sub-volume around the TERMA of each bixel
        normalizationDepth = 100 % depth [mm] of the absolute dose normalization to the machine depth dose
        kernelFile % kernel data file, default: 6 MV water kernel in matRad/doseCalc/collapsedCone
        useMex = true % use the compiled transport matRad_collapsedConeMex if available
        numThreads = 0 % number of OpenMP threads of the mex file, 0 uses the OpenMP default
        enableDijSampling = true % enable the lateral sampling of the dose influence matrix
        dijSampling % struct with the dij sampling parameters, see matRad_sampleDij
        randomSeed = 0 % seed of the random numbers used for the dij sampling
        ignoreInvalidValues = false % do not raise an error for NaN, infinite or negative doses
        kernelTilt = 'beam' % polar axis of the cones, 'beam' (beam axis) or 'ray' (ray of each bixel)
        numScatterExponentials = 1 % number of exponential terms of the scatter kernel (1 or 2)
        hardeningCorrection = true % use the fitted depth dependent effective attenuation of the TERMA
    end

    properties (SetAccess = protected, GetAccess = public)
        kernel % fitted kernel (polar bins, A, a, B, b) used for the transport
        normalizationFactor % absolute dose normalization factor
        penumbraFWHM % penumbra FWHM [mm] at the isocenter from the machine data
        mexAvailable = false % true if the compiled transport is used
        densityCubes % relative electron density cubes on the dose grid, one per CT scenario
        localDirections % cone directions in the local frame with the polar axis along z
        directionWeights % kernel weights and attenuation coefficients per component and direction
        hardeningMu0 % surface attenuation coefficient [1/mm] of the effective TERMA attenuation
        hardeningBeta = 0 % hardening coefficient [1/mm] of the effective TERMA attenuation
    end

    methods

        function this = matRad_PhotonCollapsedConeEngine(pln)
            % Constructor
            %
            % call:
            %   engine = DoseEngines.matRad_PhotonCollapsedConeEngine(pln)
            %
            % input:
            %   pln:    matRad plan meta information struct

            if nargin < 1
                pln = [];
            end

            this = this@DoseEngines.matRad_PencilBeamEngineAbstract(pln);
        end

        function setDefaults(this)
            % Sets the default properties of the engine
            setDefaults@DoseEngines.matRad_PencilBeamEngineAbstract(this);

            this.kernelFile = this.defaultKernelFile();

            this.dijSampling.relDoseThreshold = 0.01;
            this.dijSampling.latCutOff = 20;
            this.dijSampling.type = 'radius';
            this.dijSampling.deltaRadDepth = 5;
        end

    end

    methods (Access = protected)

        function dij = initDoseCalc(this, ct, cst, stf)
            % Initializes kernel, cone directions, density cubes, mex
            % availability and the absolute dose normalization
            matRad_cfg = MatRad_Config.instance();

            if any(arrayfun(@(s) ischar(s.bixelWidth), stf))
                matRad_cfg.dispError('Field based dose calculation is not supported by the %s engine.', this.name);
            end
            bixelWidth = unique([stf.bixelWidth]);
            if numel(bixelWidth) > 1
                matRad_cfg.dispError('Different bixel widths per beam are not supported!');
            end

            dij = initDoseCalc@DoseEngines.matRad_PencilBeamEngineAbstract(this, ct, cst, stf);

            if isfield(this.machine.data, 'penumbraFWHMatIso')
                this.penumbraFWHM = this.machine.data.penumbraFWHMatIso;
            else
                this.penumbraFWHM = 5;
                matRad_cfg.dispWarning('Photon machine file does not contain penumbraFWHMatIso. Assuming %g mm.', this.penumbraFWHM);
            end

            if ~any(strcmp(this.kernelTilt, {'beam', 'ray'}))
                matRad_cfg.dispError('Unknown kernel tilt ''%s'', use ''beam'' or ''ray''.', this.kernelTilt);
            end

            if this.scatterCutOff < 50
                matRad_cfg.dispWarning('Scatter cut-off of %g mm is small and truncates the scatter dose.', this.scatterCutOff);
            end

            this.effectiveLateralCutOff = max(this.geometricLateralCutOff, this.scatterCutOff) + bixelWidth / sqrt(2);

            [this.kernel, kernelMeta] = this.loadKernel(this.kernelFile, this.numPolarAngles, this.numScatterExponentials);
            [energyOk, energyMsg] = this.checkKernelEnergy(kernelMeta, this.machine);
            if ~energyOk
                matRad_cfg.dispError(energyMsg);
            elseif ~isempty(energyMsg)
                matRad_cfg.dispWarning(energyMsg);
            end
            [this.localDirections, this.directionWeights] = this.getConeDirections(this.kernel, this.numAzimuthalAngles);

            this.densityCubes = cell(1, numel(this.cubeWED));
            for s = 1:numel(this.cubeWED)
                rho = matRad_interp3(dij.ctGrid.x, dij.ctGrid.y, dij.ctGrid.z, double(this.cubeWED{s}), ...
                                     dij.doseGrid.x, dij.doseGrid.y', dij.doseGrid.z, 'linear');
                rho(~isfinite(rho) | rho < 0) = 0;
                this.densityCubes{s} = reshape(rho, dij.doseGrid.dimensions);
            end

            this.mexAvailable = this.useMex && exist('matRad_collapsedConeMex', 'file') == 3;
            if this.useMex && ~this.mexAvailable
                matRad_cfg.dispWarning(['matRad_collapsedConeMex is not compiled, using the slower MATLAB implementation. ' ...
                                        'Run matRad_compileCollapsedConeMex to compile it.']);
            end

            dz = min([dij.doseGrid.resolution.x dij.doseGrid.resolution.y dij.doseGrid.resolution.z]);
            this.hardeningMu0 = this.machine.data.m;
            this.hardeningBeta = 0;
            if this.hardeningCorrection
                [this.hardeningMu0, this.hardeningBeta] = this.fitHardening(this.kernel, this.machine, dz, this.normalizationDepth);
                matRad_cfg.dispInfo('Collapsed cone effective attenuation mu0 = %g 1/mm, beta = %g 1/mm.\n', ...
                                    this.hardeningMu0, this.hardeningBeta);
            end
            ddSvd = this.calcAnalyticalDepthDose(this.machine, this.normalizationDepth);
            ddCcc = this.calcBroadFieldDepthDose(this.kernel, this.hardeningMu0, this.normalizationDepth, dz, this.hardeningBeta);
            this.normalizationFactor = ddSvd / ddCcc;
            matRad_cfg.dispInfo('Collapsed cone normalization factor %g at %g mm depth.\n', this.normalizationFactor, this.normalizationDepth);

            if matRad_cfg.isOctave
                rand('state', this.randomSeed); %#ok<RAND>
            else
                rng(this.randomSeed);
            end
        end

        function ray = initRay(this, currBeam, j)
            % Initializes the ray geometry and the tilted cone directions in
            % dose grid array axes
            ray = initRay@DoseEngines.matRad_PencilBeamEngineAbstract(this, currBeam, j);

            rotT = double(currBeam.rotMat_system_T);
            if strcmp(this.kernelTilt, 'beam')
                rayDirWorld = [0 1 0] * rotT';
            else
                rayDirWorld = (double(ray.targetPoint_bev) - double(ray.sourcePoint_bev)) * rotT';
            end
            lateralWorld = [1 0 0] * rotT';
            ray.cccDirections = this.tiltDirections(this.localDirections, rayDirWorld, lateralWorld);
        end

        function scenRay = extractSingleScenarioRay(this, ray, scenIdx)
            % Extracts the ray of a single scenario and stores the CT
            % scenario index
            scenRay = extractSingleScenarioRay@DoseEngines.matRad_PencilBeamEngineAbstract(this, ray, scenIdx);
            scenRay.ctScen = this.multScen.linearMask(this.multScen.scenNum(scenIdx), 1);
        end

        function bixel = computeBixel(this, currRay, k)
            % Computes the collapsed cone dose of an individual bixel
            %
            % call:
            %   bixel = this.computeBixel(currRay, k)

            bixel = struct();

            if this.calcDoseDirect
                bixel.weight = currRay.weight(k);
            end

            if ~isfield(this.tmpMatrixContainers, 'physicalDose')
                bixel.ix = [];
                return
            end

            dose = this.calcBixelDose(currRay);

            if this.enableDijSampling && ~this.calcDoseDirect
                [bixel.ix, bixel.physicalDose] = matRad_sampleDij(currRay.ix, dose, currRay.radDepths, ...
                                                                  currRay.radialDist_sq, this.dijSampling, currRay.bixelWidth);
            else
                bixel.ix = currRay.ix;
                bixel.physicalDose = dose;
            end
            bixel.physicalDose = cast(bixel.physicalDose, this.precision);
        end

        function dose = calcBixelDose(this, currRay)
            % Computes TERMA, transports it on the sub-volume around the
            % bixel and returns the dose on the ray voxels
            sigma = this.penumbraFWHM / sqrt(8 * log(2));
            terma = this.calcTerma(currRay.SAD, this.hardeningMu0, sigma, currRay.bixelWidth, ...
                                   this.machine.data.primaryFluence, currRay.rayPos_bev, double(currRay.radDepths), ...
                                   double(currRay.geoDepths), double(currRay.isoLatDists), this.hardeningBeta);

            dose = zeros(numel(currRay.ix), 1);
            if ~any(terma > 0)
                return
            end

            rho = this.densityCubes{currRay.ctScen};
            dims = size(rho);
            res = [this.doseGrid.resolution.y this.doseGrid.resolution.x this.doseGrid.resolution.z];

            [i1, i2, i3] = ind2sub(dims, currRay.ix(:));
            sub = [i1 i2 i3];
            isSource = terma > 1e-6 * max(terma);
            margin = ceil(this.scatterCutOff ./ res);
            boxLo = max(min(sub(isSource, :), [], 1) - margin, 1);
            boxHi = min(max(sub(isSource, :), [], 1) + margin, dims);
            boxDims = boxHi - boxLo + 1;

            inBox = all(sub >= repmat(boxLo, size(sub, 1), 1) & sub <= repmat(boxHi, size(sub, 1), 1), 2);
            local = sub(inBox, :) - repmat(boxLo - 1, nnz(inBox), 1);
            localIx = sub2ind(boxDims, local(:, 1), local(:, 2), local(:, 3));

            termaBox = zeros(boxDims);
            termaBox(localIx) = terma(inBox);
            rhoBox = rho(boxLo(1):boxHi(1), boxLo(2):boxHi(2), boxLo(3):boxHi(3));

            w = this.directionWeights;
            doseBox = this.transport(termaBox, rhoBox, res, currRay.cccDirections, w, boxLo - 1);

            dose(inBox) = this.normalizationFactor * doseBox(localIx);

            if ~this.ignoreInvalidValues && any(~isfinite(dose) | dose < 0)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid numerical values in collapsed cone dose calculation.');
            end
        end

        function dose = transport(this, terma, rho, voxelSize, dirs, w, origin)
            % Dispatches the transport to the mex file or the MATLAB
            % implementation
            if this.mexAvailable
                dose = matRad_collapsedConeMex(terma, rho, voxelSize, dirs, w.weights, w.att, origin, this.numThreads);
            else
                dose = matRad_collapsedConeTransport(terma, rho, voxelSize, dirs, w.weights, w.att, origin);
            end
        end

    end

    methods (Static)

        function [available, msg] = isAvailable(pln, machine)
            % see superclass for information

            msg = [];
            available = false;

            if nargin < 2
                machine = matRad_loadMachine(pln);
            end

            try
                checkBasic = isfield(machine, 'meta') && isfield(machine, 'data');
                checkModality = any(strcmp(DoseEngines.matRad_PhotonCollapsedConeEngine.possibleRadiationModes, machine.meta.radiationMode));
                if ~(checkBasic && checkModality)
                    return
                end
            catch
                msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';
                return
            end

            checkData = all(isfield(machine.data, {'betas', 'm', 'primaryFluence', 'kernel', 'kernelPos'}));
            checkMeta = isfield(machine.meta, 'SAD');
            if ~(checkData && checkMeta)
                return
            end

            kernelFile = DoseEngines.matRad_PhotonCollapsedConeEngine.defaultKernelFile();
            if isstruct(pln) && isfield(pln, 'propDoseCalc') && isfield(pln.propDoseCalc, 'kernelFile')
                kernelFile = pln.propDoseCalc.kernelFile;
            end
            if exist(kernelFile, 'file') ~= 2
                msg = sprintf('Collapsed cone kernel file %s not found!', kernelFile);
                return
            end
            kernelData = load(kernelFile);
            if ~isfield(kernelData, 'ccKernel')
                msg = sprintf('Kernel file %s does not contain a collapsed cone kernel!', kernelFile);
                return
            end
            [energyOk, energyMsg] = DoseEngines.matRad_PhotonCollapsedConeEngine.checkKernelEnergy(kernelData.ccKernel.meta, machine);
            if ~energyOk
                msg = energyMsg;
                return
            end

            available = true;

            msgs = {};
            if ~isempty(energyMsg)
                msgs{end + 1} = energyMsg;
            end
            if ~isfield(machine.data, 'penumbraFWHMatIso')
                msgs{end + 1} = 'No penumbra given, generic value will be used!';
            end
            if exist('matRad_collapsedConeMex', 'file') ~= 3
                msgs{end + 1} = 'matRad_collapsedConeMex not compiled, the slower MATLAB implementation will be used!';
            end
            if ~isempty(msgs)
                msg = strjoin(msgs, ' ');
            end
        end

        function fileName = defaultKernelFile()
            % Returns the path of the default 6 MV water kernel
            fileName = fullfile(fileparts(which('matRad_fitCollapsedConeKernel')), 'matRad_collapsedConeKernel_6MV.mat');
        end

        function [ok, msg] = checkKernelEnergy(kernelMeta, machine)
            % Compares the nominal energy of the kernel with the energy of
            % the machine; a kernel without nominal energy gives a warning
            ok = true;
            msg = '';
            if ~isfield(kernelMeta, 'nominalEnergy') || ~isfield(machine.data, 'energy')
                msg = 'Nominal energy of the collapsed cone kernel or the machine unknown, energies not checked.';
                return
            end
            if abs(kernelMeta.nominalEnergy - machine.data.energy) > 1e-6
                ok = false;
                msg = sprintf(['Collapsed cone kernel for %g MV does not match the machine energy of %g MV. ' ...
                               'Generate a kernel with tools/collapsedConeKernels and set kernelFile.'], ...
                              kernelMeta.nominalEnergy, machine.data.energy);
            end
        end

        function [kernel, meta] = loadKernel(fileName, numPolarAngles, numScatterExponentials)
            % Loads the kernel data and refits it if a different number of
            % polar bins or scatter exponentials is requested
            matRad_cfg = MatRad_Config.instance();
            if nargin < 3
                numScatterExponentials = 1;
            end
            data = load(fileName, 'ccKernel');
            if ~isfield(data, 'ccKernel')
                matRad_cfg.dispError('Kernel file %s does not contain a collapsed cone kernel!', fileName);
            end
            fit = data.ccKernel.fit;
            meta = struct();
            if isfield(data.ccKernel, 'meta')
                meta = data.ccKernel.meta;
            end
            storedScatterExponentials = size(fit.B, 1);
            if fit.numPolarAngles == numPolarAngles && storedScatterExponentials == numScatterExponentials
                kernel = fit;
                kernel.numScatterExponentials = storedScatterExponentials;
                return
            end
            if ~isfield(data.ccKernel, 'cones')
                matRad_cfg.dispError('Kernel file %s cannot be refitted, it contains no conical kernel data.', fileName);
            end
            minRadius = 0;
            if isfield(fit, 'minRadius')
                minRadius = fit.minRadius;
            end
            scatterMaxRadius = Inf;
            if isfield(fit, 'scatterMaxRadius')
                scatterMaxRadius = fit.scatterMaxRadius;
            end
            kernel = matRad_fitCollapsedConeKernel(data.ccKernel.cones, numPolarAngles, minRadius, numScatterExponentials, ...
                                                   scatterMaxRadius);
        end

        function [dirs, weights] = getConeDirections(kernel, numAzimuthalAngles)
            % Returns the cone directions in the local frame (polar axis
            % along z) and the weights and attenuation coefficients of all
            % kernel components per direction (components x directions)
            [theta, phi] = ndgrid(kernel.polarAngles, ((1:numAzimuthalAngles) - 0.5) * 360 / numAzimuthalAngles);
            dirs = [sind(theta(:))' .* cosd(phi(:))'; sind(theta(:))' .* sind(phi(:))'; cosd(theta(:))'];
            weights.weights = repmat([kernel.A; kernel.B], 1, numAzimuthalAngles) / numAzimuthalAngles;
            weights.att = repmat([kernel.a; kernel.b], 1, numAzimuthalAngles);
        end

        function dirs = tiltDirections(localDirs, axisWorld, lateralWorld)
            % Rotates local cone directions so that the polar axis is along
            % axisWorld and converts them to dose grid array axes (y, x, z)
            u = axisWorld(:) / norm(axisWorld);
            e1 = lateralWorld(:) - dot(lateralWorld(:), u) * u;
            e1 = e1 / norm(e1);
            e2 = cross(u, e1);
            dirsWorld = [e1 e2 u] * localDirs;
            dirs = dirsWorld([2 1 3], :);
        end

        function terma = calcTerma(SAD, mu, sigma, bixelWidth, primaryFluence, rayPos, radDepths, geoDists, isoLatDists, beta)
            % Computes the TERMA of a bixel
            %
            % call:
            %   terma = calcTerma(SAD, mu, sigma, bixelWidth, primaryFluence, rayPos, radDepths, geoDists, isoLatDists)
            %   terma = calcTerma(SAD, mu, sigma, bixelWidth, primaryFluence, rayPos, radDepths, geoDists, isoLatDists, beta)
            %
            % input:
            %   SAD:            source to axis distance [mm]
            %   mu:             effective attenuation coefficient [1/mm]
            %   sigma:          penumbra sigma at the isocenter [mm]
            %   bixelWidth:     bixel width at the isocenter [mm]
            %   primaryFluence: radial primary fluence profile [r, value]
            %   rayPos:         ray position in beam's eye view [mm]
            %   radDepths:      radiological depths [mm]
            %   geoDists:       distances from the source [mm]
            %   isoLatDists:    lateral distances to the ray at the
            %                   isocenter plane (x, z) [mm]
            %   beta:           (optional) hardening coefficient [1/mm] of
            %                   the attenuation mu / (1 + beta d). Default: 0
            %
            % output:
            %   terma:          TERMA per unit bixel weight

            edge = @(u) 0.5 * (erf((u + bixelWidth / 2) / (sqrt(2) * sigma)) - erf((u - bixelWidth / 2) / (sqrt(2) * sigma)));
            psiLat = edge(isoLatDists(:, 1)) .* edge(isoLatDists(:, 2));

            r = sqrt((rayPos(1) + isoLatDists(:, 1)).^2 + (rayPos(3) + isoLatDists(:, 2)).^2);
            psiLat = psiLat .* interp1(primaryFluence(:, 1), primaryFluence(:, 2), r, 'linear', 0);

            if nargin < 10
                beta = 0;
            end
            terma = psiLat .* (SAD ./ geoDists(:)).^2 .* DoseEngines.matRad_PhotonCollapsedConeEngine.termaDepthFactor(mu, beta, radDepths(:));
            terma(~isfinite(terma)) = 0;
        end

        function dd = calcAnalyticalDepthDose(machine, depths, convResolution)
            % Analytical broad field central axis depth dose of the SVD
            % pencil beam machine data (without inverse square), using the
            % kernels of the SSD closest to the SAD. The kernels are
            % integrated on the same discrete grid as in the SVD pencil beam
            % engine, which resolves their central peaks.
            if nargin < 3
                convResolution = 0.5;
            end
            [~, ssdIx] = min(abs([machine.data.kernel.SSD] - machine.meta.SAD));
            kernelPos = machine.data.kernelPos(:);
            limit = ceil(kernelPos(end) / convResolution);
            [x, z] = meshgrid((-limit:limit - 1) * convResolution);
            r = sqrt(x(:).^2 + z(:).^2);
            scale = convResolution^2;
            if ~isfield(machine, 'version') || machine.version < 2
                scale = 4 * scale;
            end
            betas = machine.data.betas;
            m = machine.data.m;
            dd = zeros(size(depths));
            for i = 1:numel(betas)
                k = machine.data.kernel(ssdIx).(sprintf('kernel%d', i));
                kernelIntegral = scale * sum(interp1(kernelPos, k(:), r, 'linear', 0));
                dd = dd + kernelIntegral * betas(i) / (betas(i) - m) * (exp(-m * depths) - exp(-betas(i) * depths));
            end
        end

        function factor = termaDepthFactor(mu, beta, depths)
            % TERMA per unit primary fluence in water as a function of the
            % radiological depth for the attenuation mu / (1 + beta d)
            if beta > 0
                factor = mu ./ (1 + beta * depths) .* (1 + beta * depths).^(-mu / beta);
            else
                factor = mu * exp(-mu * depths);
            end
        end

        function [mu0, beta] = fitHardening(kernel, machine, dz, refDepth)
            % Fits mu0 and beta of the effective attenuation
            % mu0 / (1 + beta d) so that the collapsed cone broad field depth
            % dose matches the analytical SVD pencil beam depth dose of the
            % machine between 10 and 300 mm (normalized at refDepth)
            depths = (10:dz:300)';
            reference = DoseEngines.matRad_PhotonCollapsedConeEngine.calcAnalyticalDepthDose(machine, [depths; refDepth]);
            reference = reference / reference(end);
            m = machine.data.m;
            cost = @(mu, b) DoseEngines.matRad_PhotonCollapsedConeEngine.hardeningCost(kernel, mu, b, [depths; refDepth], dz, reference);
            betaGrid = [0 logspace(-5, -2, 16)];
            costGrid = zeros(size(betaGrid));
            for i = 1:numel(betaGrid)
                [~, costGrid(i)] = DoseEngines.matRad_PhotonCollapsedConeEngine.fitHardeningMu0(cost, m, betaGrid(i));
            end
            [~, iBest] = min(costGrid);
            lo = betaGrid(max(iBest - 1, 1));
            hi = betaGrid(min(iBest + 1, numel(betaGrid)));
            beta = fminbnd(@(b) DoseEngines.matRad_PhotonCollapsedConeEngine.fitHardeningMu0Cost(cost, m, b), lo, hi, optimset('TolX', 1e-7, 'Display', 'off'));
            mu0 = DoseEngines.matRad_PhotonCollapsedConeEngine.fitHardeningMu0(cost, m, beta);
        end

        function [mu0, c] = fitHardeningMu0(cost, m, beta)
            % Optimal mu0 for a given hardening coefficient
            [logMu, c] = fminbnd(@(l) cost(m * exp(l), beta), log(0.5), log(2), optimset('TolX', 1e-8, 'Display', 'off'));
            mu0 = m * exp(logMu);
        end

        function c = fitHardeningMu0Cost(cost, m, beta)
            % Cost of the optimal mu0 for a given hardening coefficient
            [~, c] = DoseEngines.matRad_PhotonCollapsedConeEngine.fitHardeningMu0(cost, m, beta);
        end

        function cost = hardeningCost(kernel, mu, beta, depths, dz, reference)
            % Relative squared deviation of the normalized broad field
            % depth dose from the reference (last depth = normalization)
            dd = DoseEngines.matRad_PhotonCollapsedConeEngine.calcBroadFieldDepthDose(kernel, mu, depths, dz, beta);
            dd = dd / dd(end);
            cost = sum(((dd(1:end - 1) - reference(1:end - 1)) ./ reference(1:end - 1)).^2);
        end

        function dd = calcBroadFieldDepthDose(kernel, mu, depths, dz, beta)
            % Collapsed cone broad field central axis depth dose in water
            % (without inverse square), computed with the recursion of the
            % transport on a one dimensional depth grid
            if nargin < 5
                beta = 0;
            end
            numVox = ceil((max(depths) + 300) / dz);
            z = ((1:numVox) - 0.5) * dz;
            terma = DoseEngines.matRad_PhotonCollapsedConeEngine.termaDepthFactor(mu, beta, z(:));
            dirs = [cosd(kernel.polarAngles); sind(kernel.polarAngles); zeros(1, kernel.numPolarAngles)];
            weights = [kernel.A; kernel.B];
            att = [kernel.a; kernel.b];
            steps = dz ./ abs(dirs(1, :));
            dose = zeros(numVox, 1);
            for forward = [true false]
                sel = (dirs(1, :) > 0) == forward;
                if ~any(sel)
                    continue
                end
                w = weights(:, sel);
                x = att(:, sel) .* repmat(steps(sel), size(att, 1), 1);
                oneMinusE = -expm1(-x);
                g = oneMinusE ./ x;
                flow = zeros(size(w));
                if forward
                    order = 1:numVox;
                else
                    order = numVox:-1:1;
                end
                for i = order
                    source = terma(i) * w;
                    dose(i) = dose(i) + sum(sum(flow .* g + source .* (1 - g)));
                    flow = flow .* (1 - oneMinusE) + source .* oneMinusE;
                end
            end
            dd = interp1(z, dose, depths, 'linear');
        end

    end
end
