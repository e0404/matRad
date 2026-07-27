classdef matRad_PhotonVmcEngine < DoseEngines.matRad_MonteCarloEngineAbstract
    % Engine for photon Monte Carlo dose calculation with VMC++.
    %   Integrates the legacy matRad_calcPhotonDoseVmc wrapper into the
    %   DoseEngines class hierarchy. The engine only reports itself as
    %   available when the (separately licensed, not distributed) VMC++
    %   environment is present under doseCalc/vmc++.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2016-2026 the matRad development team.
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
        name = 'VMC++'
        shortName = 'vmc'
    end

    properties (SetAccess = public, GetAccess = public)
        % legacy VMC++ option struct, forwarded to matRad_vmcOptions;
        % configure via pln.propDoseCalc.vmcOptions.*
        vmcOptions = struct('version', 'dkfz', ...
                            'source', 'beamlet', ...
                            'dumpDose', 2)
    end

    methods

        function this = matRad_PhotonVmcEngine(pln)
            if nargin < 1
                pln = [];
            end
            this = this@DoseEngines.matRad_MonteCarloEngineAbstract(pln);
        end

    end

    methods (Access = protected)

        function dij = calcDose(this, ct, cst, stf)
            dij = this.initDoseCalc(ct, cst, stf);

            if dij.numOfScenarios ~= 1 || ct.numOfCtScen ~= 1
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('The VMC++ dose engine currently supports only a single nominal CT scenario.');
            end

            if isfield(this.machine.data, 'weightToMU')
                dij.weightToMU = this.machine.data.weightToMU;
            else
                dij.weightToMU = 100;
            end
            dij.scaleFactor = 1;

            % assemble the legacy pln shim expected by the functional
            % VMC++ wrapper
            plnLegacy.radiationMode = stf(1).radiationMode;
            plnLegacy.propStf.numOfBeams = numel(stf);
            plnLegacy.propDoseCalc.vmcOptions = this.vmcOptions;

            dij = matRad_calcPhotonDoseVmc(ct, stf, plnLegacy, cst, this.calcDoseDirect, dij);
        end

        function dij = initDoseCalc(this, ct, cst, stf)
            % VMC++ scores on the CT voxel grid. Make that constraint
            % explicit before the common engine initializer derives grid
            % metadata and structure indices.
            ct = matRad_getWorldAxes(ct);
            requestedDoseGrid = this.doseGrid;
            useCtGrid = ~isfield(requestedDoseGrid, 'resolution') || ...
                        ~isequal(requestedDoseGrid.resolution, ct.resolution);
            if all(isfield(requestedDoseGrid, {'x', 'y', 'z'}))
                useCtGrid = useCtGrid || ~isequal(requestedDoseGrid.x, ct.x) || ...
                             ~isequal(requestedDoseGrid.y, ct.y) || ...
                             ~isequal(requestedDoseGrid.z, ct.z);
            end
            if useCtGrid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('VMC++ scores on the CT grid; the requested dose grid will be ignored.');
            end

            this.doseGrid = struct('resolution', ct.resolution, ...
                                   'x', ct.x, 'y', ct.y, 'z', ct.z);
            dij = initDoseCalc@DoseEngines.matRad_MonteCarloEngineAbstract(this, ct, cst, stf);
        end

    end

    methods (Static)

        function [available, msg] = isAvailable(pln, machine)
            % The engine requires the VMC++ environment (binaries + runs
            % folder), which is not distributed with matRad, at
            % doseCalc/vmc++ - report unavailable when it is missing.

            msg = [];
            available = false;

            if nargin < 2
                machine = matRad_loadMachine(pln);
            end

            % check basic machine file compatibility
            try
                checkBasic = isfield(machine, 'meta') && isfield(machine, 'data');
                checkModality = any(strcmp(DoseEngines.matRad_PhotonVmcEngine.possibleRadiationModes, machine.meta.radiationMode));
                if ~(checkBasic && checkModality)
                    return
                end
            catch
                msg = 'Your machine file is invalid and does not contain the basic fields (meta/data/radiationMode)!';
                return
            end

            % check for the VMC++ environment
            vmcBinDir = fullfile(DoseEngines.matRad_PhotonVmcEngine.getVmcRoot(), 'bin');
            if ~isfolder(vmcBinDir)
                msg = 'VMC++ environment not found (expected under doseCalc/vmc++/bin)!';
                return
            end

            available = true;
        end

        function vmcRoot = getVmcRoot()
            matRad_cfg = MatRad_Config.instance();
            vmcRoot = fullfile(matRad_cfg.matRadSrcRoot, 'doseCalc', 'vmc++');
        end

    end

    methods (Static, Hidden)

        function vmcCoordinates = worldToVmcCoordinates(worldCoordinates, ct)
            % The VMC CT exporter rebases the first voxel centre to one
            % voxel spacing along every axis.
            ct = matRad_getWorldAxes(ct);
            ctOriginOffset = [ct.resolution.x - ct.x(1), ...
                              ct.resolution.y - ct.y(1), ...
                              ct.resolution.z - ct.z(1)];
            vmcCoordinates = worldCoordinates + ctOriginOffset;
        end

    end
end
