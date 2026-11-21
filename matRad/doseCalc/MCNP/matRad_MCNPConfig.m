classdef matRad_MCNPConfig
    % matRad_MCNPConfig class definition
    % Holds the configurable simulation parameters of the MCNP dose engine
    %
    % References
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2024-2026 the matRad development team.
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

        %%% Simulation parameters:
        Num_Primaries = 1e6     % Number of histories per bixel
        Num_Threads             % Number of parallel calculation threads (set in constructor)
        RNG_Seed      = 43      % Seed for the random number generator

    end

    methods

        function obj = matRad_MCNPConfig()
            % matRad_MCNPConfig Configuration Class for MCNP
            matRad_cfg = MatRad_Config.instance(); % Instance of matRad configuration class

            % Number of parallel calculation threads
            if matRad_cfg.isOctave
                obj.Num_Threads = nproc();
            else
                obj.Num_Threads = feature('numcores');
            end

            % Set default histories from MatRad_Config
            if isfield(matRad_cfg.defaults.propDoseCalc, 'numHistoriesDirect')
                obj.Num_Primaries = matRad_cfg.defaults.propDoseCalc.numHistoriesDirect;
            end
        end

    end
end
