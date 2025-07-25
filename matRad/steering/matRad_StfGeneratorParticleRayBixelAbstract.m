classdef (Abstract) matRad_StfGeneratorParticleRayBixelAbstract < matRad_StfGeneratorExternalRayBixelAbstract
% matRad_StfGeneratorParticleRayBixelAbstract: Abstract Superclass for
% Particle Stf Generators using the ray-bixel mechanism
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2024 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    properties
        useRangeShifter = false;
    end

    properties (Access = protected)
        availableEnergies
        availablePeakPos
        availablePeakPosRaShi
        maxPBwidth
        pbMargin
    end

    methods
        function this = matRad_StfGeneratorParticleRayBixelAbstract(pln)
            % Constructs ExternalStfGenerator with or without pln
            if nargin < 1
                pln = [];
            end
            this@matRad_StfGeneratorExternalRayBixelAbstract(pln);

            if isempty(this.radiationMode)
                this.radiationMode = this.possibleRadiationModes{1};
            end
        end

        function setDefaults(this)
            % Set default values for ExternalStfGenerator
            this.setDefaults@matRad_StfGeneratorExternalRayBixelAbstract();
        end
    end

    methods (Access = protected)
        function initialize(this)
            this.initialize@matRad_StfGeneratorExternalRayBixelAbstract();

            %Initialize Metadata needed for stf generators
            this.availableEnergies  = [this.machine.data.energy];
            if isfield(this.machine.data, 'peakPos')
                this.availablePeakPos   = [this.machine.data.peakPos] + [this.machine.data.offset];
            else
                this.availablePeakPos   = [];
            end
            availableWidths         = [this.machine.data.initFocus];
            availableWidths         = [availableWidths.SisFWHMAtIso];
            this.maxPBwidth         = max(availableWidths) / 2.355;

            matRad_cfg = MatRad_Config.instance();
            if this.useRangeShifter
                %For now only a generic range shifter is used whose thickness is
                %determined by the minimum peak width to play with
                rangeShifterEqD = round(min(this.availablePeakPos)* 1.25);
                this.availablePeakPosRaShi = this.availablePeakPos - rangeShifterEqD;

                matRad_cfg.dispWarning('Use of range shifter enabled. matRad will generate a generic range shifter with WEPL %f to enable ranges below the shortest base data entry.',rangeShifterEqD);
            end

            if sum(this.availablePeakPos<0)>0
                matRad_cfg.dispError('at least one available peak position is negative - inconsistent machine file')
            end

            %Create Water equivalent cube in ct
            this.ct = matRad_calcWaterEqD(this.ct,this.radiationMode);
        end
    end

    methods (Static)
        function [available,msg] = isAvailable(pln,machine)
            % see superclass for information

            if nargin < 2
                machine = matRad_loadMachine(pln);
            end

            % Check superclass availability
            [available,msg] = matRad_StfGeneratorExternalRayBixelAbstract.isAvailable(pln,machine);

            if ~available
                return;
            end

            available = available && isstruct(machine.data);
            available = available && all(isfield(machine.data,{'energy','initFocus'}));

            if ~available
                msg = 'Your machine file is invalid and does not contain the basic fields required for photon machines!';
            else
                msg = [];
            end
        end
    end
end

