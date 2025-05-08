classdef matRad_MCNPConfig
% matRad_MCNPConfig class definition
% 
%
% References
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2024
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    properties              

        %%% Simulation parameters:
        Num_Primaries = 1e6;
        Num_Threads   =	feature('numcores');		% Number of parallel calculation threads
        RNG_Seed      =	43;		% Seed for the random number generator

        

        
    end
    
    methods
        function obj = matRad_MCNPConfig()
            %matRad_MCNPConfig Configuration Class for MCNP   
            matRad_cfg = MatRad_Config.instance(); %Instance of matRad configuration class
            
            % Set default histories from MatRad_Config
            if isfield(matRad_cfg.propDoseCalc,'defaultNumHistories')
                obj.Num_Primaries = matRad_cfg.propMC.defaultNumHistories;
            end
        end
    end
end

