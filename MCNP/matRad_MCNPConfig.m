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
        Num_Primaries = 1e6;

        %%% Simulation parameters:
        Num_Threads   =	feature('numcores');		% Number of parallel calculation threads. Default: 0 = max available threads
        RNG_Seed      =	43;		% Seed for the random number generator
        
        % % This parameter can be overwritten through MatRad_Config default parameters
        % E_Cut_Pro     =	0.5;		% Energy cut (in MeV) below which heavy charged particles are locally absorbed. Default: 0.5
    end
    
    methods
        function obj = matRad_MCNPConfig()
            %matRad_MCNPConfig Configuration Class for MCNP   
            matRad_cfg = MatRad_Config.instance(); %Instance of matRad configuration class
            
            % Set default histories from MatRad_Config
            if isfield(matRad_cfg.propMC,'defaultNumHistories')
                obj.Num_Primaries = matRad_cfg.propMC.defaultNumHistories;
            end
        end
    end
end

