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
    %     numPolarAngles      number of polar cone bins (divisor of 48)
    %     numAzimuthalAngles  number of azimuthal cones per polar bin
    %     scatterCutOff       lateral margin [mm] of the transport sub-volume
    %                         around the TERMA of each bixel. Dose is
    %                         stored within max(geometricLateralCutOff,
    %                         scatterCutOff) around each ray
    %     normalizationDepth  depth [mm] of the absolute dose normalization
    %     kernelFile          kernel data file (default: 6 MV water kernel)
    %     useMex              use the compiled transport if available
    %     numThreads          number of OpenMP threads of the mex file
    %                         (0 = default)
    %     kernelTilt          polar axis of the cones: 'beam' (beam axis,
    %                         same transport lattice for all bixels of a
    %                         beam) or 'ray' (ray of each bixel, which
    %                         causes lattice aliasing in the summed dose)
    %     enableDijSampling   enable the lateral dij sampling
    %     dijSampling         sampling parameters, see matRad_sampleDij
    %
    %   Limitations: the TERMA uses a single effective attenuation
    %   coefficient (no spectral hardening or off-axis softening) and the
    %   kernel is computed in water and scaled with the relative electron
    %   density.
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
        possibleRadiationModes = {'photons'}
        name = 'Collapsed Cone Convolution'
        shortName = 'CCC'
    end

    properties (SetAccess = public, GetAccess = public)
        numPolarAngles = 12
        numAzimuthalAngles = 8
        scatterCutOff = 150
        normalizationDepth = 100
        kernelFile
        useMex = true
        numThreads = 0
        enableDijSampling = true
        dijSampling
        randomSeed = 0
        ignoreInvalidValues = false
        kernelTilt = 'beam'
    end

    properties (SetAccess = protected, GetAccess = public)
        kernel
        normalizationFactor
        penumbraFWHM
        mexAvailable = false
        densityCubes
        localDirections
        directionWeights
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

            this.kernel = this.loadKernel(this.kernelFile, this.numPolarAngles);
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

            ddSvd = this.calcAnalyticalDepthDose(this.machine, this.normalizationDepth);
            ddCcc = this.calcBroadFieldDepthDose(this.kernel, this.machine.data.m, this.normalizationDepth, ...
                                                 min([dij.doseGrid.resolution.x dij.doseGrid.resolution.y dij.doseGrid.resolution.z]));
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
            terma = this.calcTerma(currRay.SAD, this.machine.data.m, sigma, currRay.bixelWidth, ...
                                   this.machine.data.primaryFluence, currRay.rayPos_bev, double(currRay.radDepths), ...
                                   double(currRay.geoDepths), double(currRay.isoLatDists));

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
                dose = matRad_collapsedConeMex(terma, rho, voxelSize, dirs, w.A, w.a, w.B, w.b, origin, this.numThreads);
            else
                dose = matRad_collapsedConeTransport(terma, rho, voxelSize, dirs, w.A, w.a, w.B, w.b, origin);
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

            available = true;

            msgs = {};
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

        function kernel = loadKernel(fileName, numPolarAngles)
            % Loads the kernel data and refits it if a different number of
            % polar bins is requested
            matRad_cfg = MatRad_Config.instance();
            data = load(fileName, 'ccKernel');
            if ~isfield(data, 'ccKernel')
                matRad_cfg.dispError('Kernel file %s does not contain a collapsed cone kernel!', fileName);
            end
            if data.ccKernel.fit.numPolarAngles == numPolarAngles
                kernel = data.ccKernel.fit;
            else
                minRadius = 0;
                if isfield(data.ccKernel.fit, 'minRadius')
                    minRadius = data.ccKernel.fit.minRadius;
                end
                kernel = matRad_fitCollapsedConeKernel(data.ccKernel.cones, numPolarAngles, minRadius);
            end
        end

        function [dirs, weights] = getConeDirections(kernel, numAzimuthalAngles)
            % Returns the cone directions in the local frame (polar axis
            % along z) and the kernel parameters per direction
            [theta, phi] = ndgrid(kernel.polarAngles, ((1:numAzimuthalAngles) - 0.5) * 360 / numAzimuthalAngles);
            dirs = [sind(theta(:))' .* cosd(phi(:))'; sind(theta(:))' .* sind(phi(:))'; cosd(theta(:))'];
            perDir = @(v) reshape(repmat(v(:), 1, numAzimuthalAngles), 1, []);
            weights.A = perDir(kernel.A) / numAzimuthalAngles;
            weights.a = perDir(kernel.a);
            weights.B = perDir(kernel.B) / numAzimuthalAngles;
            weights.b = perDir(kernel.b);
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

        function terma = calcTerma(SAD, mu, sigma, bixelWidth, primaryFluence, rayPos, radDepths, geoDists, isoLatDists)
            % Computes the TERMA of a bixel
            %
            % call:
            %   terma = calcTerma(SAD, mu, sigma, bixelWidth, primaryFluence, rayPos, radDepths, geoDists, isoLatDists)
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
            %
            % output:
            %   terma:          TERMA per unit bixel weight

            edge = @(u) 0.5 * (erf((u + bixelWidth / 2) / (sqrt(2) * sigma)) - erf((u - bixelWidth / 2) / (sqrt(2) * sigma)));
            psiLat = edge(isoLatDists(:, 1)) .* edge(isoLatDists(:, 2));

            r = sqrt((rayPos(1) + isoLatDists(:, 1)).^2 + (rayPos(3) + isoLatDists(:, 2)).^2);
            psiLat = psiLat .* interp1(primaryFluence(:, 1), primaryFluence(:, 2), r, 'linear', 0);

            terma = psiLat .* (SAD ./ geoDists(:)).^2 .* mu .* exp(-mu * radDepths(:));
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

        function dd = calcBroadFieldDepthDose(kernel, mu, depths, dz)
            % Collapsed cone broad field central axis depth dose in water
            % (without inverse square), computed with the same transport on
            % a one dimensional depth grid
            numVox = ceil((max(depths) + 300) / dz);
            z = ((1:numVox) - 0.5) * dz;
            terma = mu * exp(-mu * z(:));
            dirs = [cosd(kernel.polarAngles); sind(kernel.polarAngles); zeros(1, kernel.numPolarAngles)];
            dose = matRad_collapsedConeTransport(terma, ones(numVox, 1), [dz 1e9 1e9], dirs, ...
                                                 kernel.A, kernel.a, kernel.B, kernel.b);
            dd = interp1(z, dose, depths, 'linear');
        end

    end
end
