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
        % the actual calculation method wich returns the final dij struct.
        % Needs to be implemented in non abstract subclasses. 
        %(Internal logic is often split into multiple methods in order to make the whole calculation more modular)
        function dij = calcDose(this,ct,stf,pln,cst)
            error('Function needs to be implemented!');
        end
        
    end 
    
    methods(Static)
        [nameList, classList, handleList] = getAvailableEngines(pln,optionalPath)         

        function [available,msg] = isAvailable(pln,machine)   
        % return a boolean if the engine is is available for the given pln
        % struct. Needs to be implemented in non abstract subclasses
        % input:
        % - pln:        matRad pln struct
        % - machine:    optional machine to avoid loading the machine from
        %               disk (makes sense to use if machine already loaded)
        % output:
        % - available:  boolean value to check if the dose engine is 
        %               available for the given pln/machine
        % - msg:        msg to elaborate on availability. If not available,
        %               a msg string indicates an error during the check
        %               if available, indicates a warning that not all
        %               information was present in the machine file and
        %               approximations need to be made
        
            error('This is an Abstract Base class! Function needs to be called for instantiable subclasses!');
        end
    end


end
