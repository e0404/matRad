classdef (Abstract) matRad_MonteCarloEngine < DoseCalcEngines.matRad_DoseCalcEngine
    % Superclass for all dose calculation engines which are based on 
    % monte carlo calculation 
    % for more informations see superclass
    % DoseCalcEngines.matRad_DoseCalcEngine
    
    properties (SetAccess = protected, GetAccess = public)
        
        pbCalcMode; % fine sampling mode
  
        geoDistVdoseGrid;   % geometric distance in dose grid
        rot_coordsVdoseGrid;    % Rotate coordinates for gantry movement
        radDepthsMat;   % radiological depth cube container
        radDepthVdoseGrid;  % grid for radiologica depth cube
       
        effectiveLateralCutoff; % lateral cutoff for raytracing and geo calculations
        
        bixelsPerBeam;  % number of bixel per energy beam

        nCasePerBixel; %number of histories per beamlet
        
    end
    
end

