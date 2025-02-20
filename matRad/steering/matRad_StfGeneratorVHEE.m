classdef matRad_StfGeneratorVHEE < matRad_StfGeneratorParticleRayBixelAbstract
    % matRad_StfGeneratorVHEE: STF generator for Very High Energy Electrons (VHEE)
    %
    % Implements a simplified single-energy, single-spot-per-ray approach, 
    % ignoring multi-layer logic typical for scanned protons/ions.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2024 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No
    % part of the matRad project, including this file, may be copied,
    % modified, propagated, or distributed except according to the terms
    % contained in the LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Constant)
        name                   = 'VHEE STF Generator'
        shortName              = 'ParticleVHEE'
        possibleRadiationModes = {'vhee'}
    end

    %% ----------------------------------------------------------------
    %  (1) Public Methods: Constructor, setDefaults, etc.
    %% ----------------------------------------------------------------
    methods
        function this = matRad_StfGeneratorVHEE(pln)
            % Constructor
            if nargin < 1
                pln = [];
            end
            % Call the ParticleRayBixelAbstract constructor
            this@matRad_StfGeneratorParticleRayBixelAbstract(pln);

            % If plan has no radiationMode, set to 'vhee'
            if isempty(this.radiationMode)
                this.radiationMode = 'vhee';
            end
        end

        function setDefaults(this)
            % Must be public if parent is public
            setDefaults@matRad_StfGeneratorParticleRayBixelAbstract(this);
            % Additional VHEE defaults could go here
        end
    end

    %% ----------------------------------------------------------------
    %  (2) Protected Methods: Overriding initialize, plus Beam Data, etc.
    %% ----------------------------------------------------------------
    methods (Access = protected)
        
        % --- FIX A: Overriding the parent's initialize() ---
        function initialize(this)
            % Here we SKIP the ParticleRayBixelAbstract.initialize, which references peakPos.
            % Instead, we only call the next-level ExternalRayBixelAbstract.initialize
            % (or you could skip calling it entirely if your VHEE generator doesn't need it).

            this.initialize@matRad_StfGeneratorExternalRayBixelAbstract();

            % Add any additional checks or initializations for VHEE here, e.g.:
            % matRad_cfg = MatRad_Config.instance();
            % if ~isfield(this.machine.data,'energy') && ~isfield(this.machine.data,'energies')
            %     matRad_cfg.dispError('VHEE machine must have "energy" or "energies"!');
            % end
        end

        function beam = initBeamData(this, beam)
            % Inherit standard logic from the parent
            beam = initBeamData@matRad_StfGeneratorParticleRayBixelAbstract(this, beam);

            % If user doesn't specify an energy, default to 200 MeV
            if ~isfield(this.pln,'energy') || isempty(this.pln.energy)
                beam.VHEEenergy = 200; 
            else
                beam.VHEEenergy = this.pln.energy;
            end

            % If 'machine.data.energies' is present, verify it
            if isfield(this.machine.data,'energies') && ~isempty(this.machine.data.energies)
                if ~ismember(beam.VHEEenergy, this.machine.data.energies)
                    error(['The specified VHEE energy (',num2str(beam.VHEEenergy), ...
                           ' MeV) is not in machine.data.energies!']);
                end
            else
                warning(['Machine base data does not contain "energies". Using user-defined energy: ', ...
                         num2str(beam.VHEEenergy),' MeV.']);
            end
        end

        function beam = initRays(this, beam)
            beam = initRays@matRad_StfGeneratorParticleRayBixelAbstract(this, beam);
            % Additional geometry modifications, if any
        end

        function beam = setBeamletEnergies(this, beam)
            % Assign the single VHEEenergy to each ray
            for iRay = 1:numel(beam.ray)
                beam.ray(iRay).energy = beam.VHEEenergy;
            end
        end

        function beam = createRayBixels(this, beam)
            % Create exactly one bixel per ray
            for iRay = 1:numel(beam.ray)
                bixel.energy = beam.VHEEenergy;
                bixel.weight = 1.0; 
                bixel.idx    = iRay; 
                beam.ray(iRay).bixel = bixel;
            end
        end
    end

    %% ----------------------------------------------------------------
    %  (3) Overridden isAvailable: Skip Proton-Specific Checks
    %% ----------------------------------------------------------------
    methods (Static)
        function [available,msg] = isAvailable(pln,machine)
            % We do NOT call ParticleRayBixelAbstract.isAvailable
            % but jump to ExternalRayBixelAbstract to skip 'peakPos' etc.

            if nargin < 2
                machine = matRad_loadMachine(pln);
            end

            [available,msg] = matRad_StfGeneratorExternalRayBixelAbstract.isAvailable(pln,machine);
            if ~available
                return;
            end

            if ~ismember(pln.radiationMode, matRad_StfGeneratorVHEE.possibleRadiationModes)
                available = false;
                msg = ['This generator only supports radiationMode = ''vhee'', not ''', pln.radiationMode, '''.'];
                return;
            end

            if ~isfield(machine,'data') || ...
               (~isfield(machine.data,'energy') && ~isfield(machine.data,'energies'))
                available = false;
                msg = 'Machine data must contain "energy" or "energies" for VHEE.';
                return;
            end

            msg = [];
        end
    end
end