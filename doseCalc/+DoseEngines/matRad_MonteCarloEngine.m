classdef (Abstract) matRad_MonteCarloEngine < DoseEngines.matRad_DoseEngine
    % Superclass for all dose calculation engines which are based on 
    % monte carlo calculation 
    % for more informations see superclass
    % DoseEngines.matRad_DoseEngine
    
    properties (SetAccess = public, GetAccess = public)
                         
        nCasePerBixel; %number of histories per beamlet
        
    end
    
    properties (Constant)
        
        isMCEngine = true; %constant to differentiate pencil beam and mc engines
        
    end
    
    methods(Access = protected)
         
    end
    
end

