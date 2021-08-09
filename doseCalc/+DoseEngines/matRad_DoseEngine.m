classdef (Abstract) matRad_DoseEngine < handle
    % matRad_DoseEngine: Interface for dose calculation
    %   This base class provides the structure for the basic initialization 
    %   functions and corresponding properties for e.g. particle and photon
    %   based dose calc.
    %   Implementations like particle, photon based dose calculation can be
    %   found in the DoseEngine package
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        calcDoseDirect = false; % analytical mode
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
        function obj = matRad_DoseEngine()
        % future code for property validation on creation here
        end
    end
    
    methods(Access  =  protected)
        
        function [ct,stf,pln,dij] = calcDoseInit(obj,ct,stf,pln,cst)
        % method for setting and preparing the inition parameters for the 
        % dose calculation.
        % Shoulb be called at the beginning of calcDose method.
        % Can be expanded or changed by overwriting this method and calling
        % the superclass method inside of it
       
            matRad_cfg =  MatRad_Config.instance();

            % to guarantee downwards compatibility with data that does not have
            % ct.x/y/z
            if ~any(isfield(ct,{'x','y','z'}))
                ct.x = ct.resolution.x*[0:ct.cubeDim(2)-1]-ct.resolution.x/2;
                ct.y = ct.resolution.y*[0:ct.cubeDim(1)-1]-ct.resolution.y/2;
                ct.z = ct.resolution.z*[0:ct.cubeDim(3)-1]-ct.resolution.z/2;
            end

            % set grids
            if ~isfield(pln,'propDoseCalc') || ...
               ~isfield(pln.propDoseCalc,'doseGrid') || ...
               ~isfield(pln.propDoseCalc.doseGrid,'resolution')
                % default values
                dij.doseGrid.resolution = matRad_cfg.propDoseCalc.defaultResolution;
            else
                % take values from pln strcut
                dij.doseGrid.resolution.x = pln.propDoseCalc.doseGrid.resolution.x;
                dij.doseGrid.resolution.y = pln.propDoseCalc.doseGrid.resolution.y;
                dij.doseGrid.resolution.z = pln.propDoseCalc.doseGrid.resolution.z;
            end

            dij.doseGrid.x = ct.x(1):dij.doseGrid.resolution.x:ct.x(end);
            dij.doseGrid.y = ct.y(1):dij.doseGrid.resolution.y:ct.y(end);
            dij.doseGrid.z = ct.z(1):dij.doseGrid.resolution.z:ct.z(end);

            dij.doseGrid.dimensions  = [numel(dij.doseGrid.y) numel(dij.doseGrid.x) numel(dij.doseGrid.z)];
            dij.doseGrid.numOfVoxels = prod(dij.doseGrid.dimensions);

            dij.ctGrid.resolution.x = ct.resolution.x;
            dij.ctGrid.resolution.y = ct.resolution.y;
            dij.ctGrid.resolution.z = ct.resolution.z;

            dij.ctGrid.x = ct.x;
            dij.ctGrid.y = ct.y;
            dij.ctGrid.z = ct.z;

            dij.ctGrid.dimensions  = [numel(dij.ctGrid.y) numel(dij.ctGrid.x) numel(dij.ctGrid.z)];
            dij.ctGrid.numOfVoxels = prod(dij.ctGrid.dimensions);

            % adjust isocenter internally for different dose grid
            obj.offset = [dij.doseGrid.resolution.x - dij.ctGrid.resolution.x ...
                      dij.doseGrid.resolution.y - dij.ctGrid.resolution.y ...
                      dij.doseGrid.resolution.z - dij.ctGrid.resolution.z];

            for i = 1:numel(stf)
                stf(i).isoCenter = stf(i).isoCenter + obj.offset;
            end

            %set up HU to rED or rSP conversion
            if ~isfield(pln,'propDoseCalc') || ~isfield(pln.propDoseCalc,'useGivenEqDensityCube')
                disableHUconversion = matRad_cfg.propDoseCalc.defaultUseGivenEqDensityCube;
            else
                disableHUconversion = pln.propDoseCalc.useGivenEqDensityCube;
            end

            %If we want to omit HU conversion check if we have a ct.cube ready
            if disableHUconversion && ~isfield(ct,'cube')
                matRad_cfg.dispWarning('HU Conversion requested to be omitted but no ct.cube exists! Will override and do the conversion anyway!');
                disableHUconversion = false;
            end

            % calculate rED or rSP from HU
            if disableHUconversion
                matRad_cfg.dispInfo('Omitting HU to rED/rSP conversion and using existing ct.cube!\n');
            else
                ct = matRad_calcWaterEqD(ct, pln);
            end

            % meta information for dij
            dij.numOfBeams         = pln.propStf.numOfBeams;
            dij.numOfScenarios     = 1;
            dij.numOfRaysPerBeam   = [stf(:).numOfRays];
            dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
            dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam); 
            
            % check if full dose influence data is required
            if obj.calcDoseDirect 
                obj.numOfColumnsDij      = length(stf);
                obj.numOfBixelsContainer = 1;
            else
                obj.numOfColumnsDij      = dij.totalNumOfBixels;
                obj.numOfBixelsContainer = ceil(dij.totalNumOfBixels/10);
            end

            % set up arrays for book keeping
            dij.bixelNum = NaN*ones(obj.numOfColumnsDij,1);
            dij.rayNum   = NaN*ones(obj.numOfColumnsDij,1);
            dij.beamNum  = NaN*ones(obj.numOfColumnsDij,1);


            % Allocate space for dij.physicalDose sparse matrix
            for i = 1:dij.numOfScenarios
                dij.physicalDose{i} = spalloc(dij.doseGrid.numOfVoxels,obj.numOfColumnsDij,1);
            end

            % take only voxels inside patient
            VctGrid = [cst{:,4}];
            VctGrid = unique(vertcat(VctGrid{:}));

            % ignore densities outside of contours
            if ~isfield(pln,'propDoseCalc') || ~isfield(pln.propDoseCalc,'ignoreOutsideDensities')
                ignoreOutsideDensities = matRad_cfg.propDoseCalc.defaultIgnoreOutsideDensities;
            else
                ignoreOutsideDensities = pln.propDoseCalc.ignoreOutsideDensities;
            end

            if ignoreOutsideDensities
                eraseCtDensMask = ones(prod(ct.cubeDim),1);
                eraseCtDensMask(VctGrid) = 0;
                for i = 1:ct.numOfCtScen
                    ct.cube{i}(eraseCtDensMask == 1) = 0;
                end
            end                       
            
            % receive linear indices and grid locations from the dose grid
            tmpCube    = zeros(ct.cubeDim);
            tmpCube(VctGrid) = 1;
            % interpolate cube
            obj.VdoseGrid = find(matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y,   dij.ctGrid.z,tmpCube, ...
                                            dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'nearest'));
                                                                              
            % save vct grid as own property in order to allow sub-classes
            % to access it
            obj.VctGrid = VctGrid;
            
            % Convert CT subscripts to linear indices.
            [obj.yCoordsV_vox, obj.xCoordsV_vox, obj.zCoordsV_vox] = ind2sub(ct.cubeDim,obj.VctGrid);
            
            
            % Convert CT subscripts to coarse linear indices.
            [obj.yCoordsV_voxDoseGrid, obj.xCoordsV_voxDoseGrid, obj.zCoordsV_voxDoseGrid] = ind2sub(dij.doseGrid.dimensions,obj.VdoseGrid);

            % load machine file from base data folder
            obj.machine = obj.loadMachine(pln);
            
            % compute SSDs
            stf = matRad_computeSSD(stf,ct);

        end
                
    end
    
    methods
        
        function set.calcDoseDirect(obj,calcDoseDirect)
%             %SETTER with argument validation for calcDoseDirect property
%             Sadly we can't use the argument validation of matlab because
%             octave doesn't support this
%             arguments     
%                 obj
%                 calcDoseDirect {mustBeOfClass(calcDoseDirect, 'logical')}           
%             end  
%            
            
            obj.calcDoseDirect = calcDoseDirect;
        end
        
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
        function dij = calcDose(obj,ct,stf,pln,cst)
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
        
    end
end
