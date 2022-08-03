classdef (Abstract) matRad_DoseEngine < handle
% matRad_DoseEngine: Interface for dose calculation
%   This base class provides the structure for the basic initialization 
%   functions and corresponding properties for e.g. particle and photon
%   based dose calc.
%   Implementations like particle, photon based dose calculation can be
%   found in the DoseEngine package
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        machine; % base data defined in machine file
    end
    
    properties (SetAccess = protected, GetAccess = public)
        
        numOfBixelsContainer;   % number of used bixel container
        numOfColumnsDij;    % number of columns in the dij struct
                                          
        yCoordsV_vox;   % y-coordinate voxel
        xCoordsV_vox;   % x-coordinate voxel
        zCoordsV_vox;   % z-coordinate voxel
        
        yCoordsV_voxDoseGrid;   % converted voxel indices to real grid 
        xCoordsV_voxDoseGrid;   % converted voxel indices to real grid
        zCoordsV_voxDoseGrid;   % converted voxel indices to real grid
        
        offset; % offset adjustment for isocenter
        
        VctGrid; % voxel grid inside patient
        VdoseGrid;  % voxel dose grid 
        
    end
  
    properties (SetAccess = public, GetAccess = public)
        calcDoseDirect = false; % switch for direct cube / dij calculation 
    end
    
    properties (Constant)
        isDoseEngine = true; % const boolean for checking inheritance
    end
    
    properties (Constant, Abstract)
       name; % readable name for dose engine
       possibleRadiationModes; % radiation modes the engine is meant to process
    end
    
    methods      
        %Constructor  
        function this = matRad_DoseEngine()
        % future code for property validation on creation here
        end
    end
    
    methods(Access  =  protected)
        
        % method for setting and preparing the inition parameters for the 
        % dose calculation.
        % Should be called at the beginning of calcDose method.
        % Can be expanded or changed by overwriting this method and calling
        % the superclass method inside of it
        [ct,stf,pln,dij] = calcDoseInit(this,ct,cst,pln,stf)
        
        
    end
    
    % Should be abstract methods but in order to satisfy the compatibility
    % with OCTAVE we can't use abstract methods. If OCTAVE at some point 
    % in the far future implements this feature this should be abstract again.
    methods %(Abstract)
        
       
        function ret = isAvailable(pln)   
        % return a boolean if the engine is is available for the given pln
        % struct. Needs to be implemented in non abstract subclasses
            error('Funktion needs to be implemented!');
        end
        
        % the actual calculation method wich returns the final dij struct.
        % Needs to be implemented in non abstract subclasses. 
        %(Internal logic is often split into multiple methods in order to make the whole calculation more modular)
        function dij = calcDose(this,ct,stf,pln,cst)
            error('Funktion needs to be implemented!');
        end
        
    end 
    
    methods(Static)
        
        function machine = loadMachine(pln, filepath)
            %load the machine mode for the specific dose calculation
            %static so it can be used outside e.g. to validate the machine,
            %pre the dose calculation
            matRad_cfg = MatRad_Config.instance();
            if isfield(pln, 'radiationMode')
                if isfield(pln, 'machine')
                    fileName = [pln.radiationMode '_' pln.machine];
                else
                    fileName = [pln.radiationMode '_Generic'];
                end 
            else
                matRad_cfg.dispError('No radiation mode given in pln');
            end 
            
            if ~exist('filepath','var')
                filepath = [matRad_cfg.matRadRoot filesep 'basedata' filesep fileName];
            end
            
            try
               m = load(filepath, 'machine');
               machine = m.machine; % strip first layer of loaded struct for convenience 
            catch
               matRad_cfg.dispWarning('Could not find the following machine file: %s\n',fileName); 
            end
        end

        [nameList, classList, handleList] = getAvailableEngines(pln,optionalPath)     

        [nameList, classList, handleList] = getAvailableEnginesOctave(pln,optionalPath)     
    end
end
