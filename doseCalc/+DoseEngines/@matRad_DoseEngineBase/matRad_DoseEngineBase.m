classdef (Abstract) matRad_DoseEngineBase < handle
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
    
    % Implement these properties in a subclass
    properties (Constant, Abstract)
       shortName;               % short identifier by which matRad recognizes an engine
       name;                    % user readable name for dose engine
       possibleRadiationModes;  % radiation modes the engine is meant to process
       %supportedQuantities;    % supported (influence) quantities. Does not include quantities that can be derived post-calculation.
    end    
    
    % Public properties
    properties
        doseGrid;                   % doseGrid to use (struct with at least doseGrid.resolution.x/y/z set)        
        multScen;                   % scenario model to use
        voxelSubIx;                 % selection of where to calculate / store dose, empty by default
        selectVoxelsInScenarios;    % which voxels to compute in robustness scenarios
        bioParam;                   % Biological dose modeling
    end
    
    % Protected properties with public get access
    properties (SetAccess = protected, GetAccess = public)
        machine;                % base data defined in machine file

        timers;                 % timers of dose calc

        numOfColumnsDij;        % number of columns in the dij struct
                                          
        yCoordsV_vox;           % y-coordinate voxel
        xCoordsV_vox;           % x-coordinate voxel
        zCoordsV_vox;           % z-coordinate voxel
        
        yCoordsV_voxDoseGrid;   % converted voxel indices to real grid 
        xCoordsV_voxDoseGrid;   % converted voxel indices to real grid
        zCoordsV_voxDoseGrid;   % converted voxel indices to real grid
        
        %offset; % offset adjustment for isocenter
        
        VctGrid; % voxel grid inside patient
        VdoseGrid;  % voxel dose grid    
        VctGridScenIx;  %logical subindexes of scenarios in ct grid
        VdoseGridScenIx; %logical subindexes of scenarios in dose grid

        VctGridMask; % voxel grid inside patient as logical mask
        VdoseGridMask;  % voxel dose grid inside patient as logical mask

        robustVoxelsOnGrid; %voxels to be computed in robustness scenarios
    end
    
    % Fully protected properties
    properties (Access = protected)
        lastProgressUpdate;
        calcDoseDirect = false; % switch for direct cube / dij calculation
    end
    
    properties (Constant)
        isDoseEngine = true;    % const boolean for checking inheritance
    end
    
    properties (SetAccess = private)
        hWaitbar;
    end
    
    methods      
        %Constructor  
        function this = matRad_DoseEngineBase(pln)
            this.setDefaults();
            if nargin == 1 && ~isempty(pln)
                this.assignPropertiesFromPln(pln);
            end
        end

        function warnDeprecatedEngineProperty(this,oldProp,msg,newProp)
            matRad_cfg = MatRad_Config.instance();
            if nargin < 3 || isempty(msg)
                msg = '';
            end

            if nargin < 4
                dep2 = '';
            else
                dep2 = sprintf('Use Property ''%s'' instead!',newProp);
            end

            matRad_cfg.dispDeprecationWarning('Property ''%s'' of Dose Engine ''%s'' is deprecated! %s%s',oldProp,this.name,msg,dep2);
        end

        function assignPropertiesFromPln(this,pln,warnWhenPropertyChanged)
            matRad_cfg = MatRad_Config.instance();
            
            %Set Scenario Model
            if isfield(pln,'multScen')
                this.multScen = pln.multScen;
            end
            
            %Assign biological model
            if isfield(pln,'bioParam')
                this.bioParam = pln.bioParam;
            end
            
            if nargin < 3 || ~isscalar(warnWhenPropertyChanged) || ~islogical(warnWhenPropertyChanged)
                warnWhenPropertyChanged = false;
            end

            %Overwrite default properties within the engine with the ones
            %given in the propDoseCalc struct
            if isfield(pln,'propDoseCalc') && isstruct(pln.propDoseCalc)
                fields = fieldnames(pln.propDoseCalc); %get remaining fields
                if isfield(pln.propDoseCalc,'engine') && ~isempty(pln.propDoseCalc.engine) && ~strcmp(pln.propDoseCalc.engine,this.shortName)
                    matRad_cfg.dispWarning('Inconsistent dose engines given! pln asks for ''%s'', but you are using ''%s''!',pln.propDoseCalc.engine,this.shortName);
                end
                fields(strcmp(fields, 'engine')) = []; % engine field is no longer needed and would throw an exception
            else
                fields = {};
            end

            % iterate over all fieldnames and try to set the
            % corresponding properties inside the engine
            for i = 1:length(fields)
                try
                    oldValue = this.(fields{i});
                    newValue = pln.propDoseCalc.(fields{i});
                    this.(fields{i}) = newValue;

                    if warnWhenPropertyChanged
                        if ~isequal(oldValue,newValue)
                            matRad_cfg.dispWarning('Property ''%s'' overwritten by Plan settings!',fields{i});
                        end
                    end

                % catch exceptions when the engine has no properties,
                % which are defined in the struct.
                % When defining an engine with custom setter and getter
                % methods, custom exceptions can be caught here. Be
                % careful with Octave exceptions!
                catch ME
                    switch ME.identifier
                        case 'MATLAB:noPublicFieldForClass'
                            matRad_cfg.dispWarning('Not able to assign property from pln.propDoseCalc to Dose Engine: %s',ME.message);
                        otherwise
                            matRad_cfg.dispWarning('Problem while setting up engine from struct:%s %s',fields{i},ME.message);
                    end
                end
            end           
        end
    
        function resultGUI = calcDoseForward(this,ct,cst,stf,w)
            matRad_cfg = MatRad_Config.instance();
            if nargin < 5 && ~isfield([stf.ray],'weight')
                matRad_cfg.dispError('No weight vector available. Please provide w or add info to stf')
            end

            % copy bixel weight vector into stf struct
            if nargin == 5
                if sum([stf.totalNumOfBixels]) ~= numel(w)
                    matRad_cfg.dispError('weighting does not match steering information');
                end
                counter = 0;
                for i = 1:size(stf,2)
                    for j = 1:stf(i).numOfRays
                        for k = 1:stf(i).numOfBixelsPerRay(j)
                            counter = counter + 1;
                            stf(i).ray(j).weight(k) = w(counter);
                        end
                    end
                end
            else % weights need to be in stf!
                w = NaN*ones(sum([stf.totalNumOfBixels]),1);
                counter = 0;
                for i = 1:size(stf,2)
                    for j = 1:stf(i).numOfRays
                        for k = 1:stf(i).numOfBixelsPerRay(j)
                            counter = counter + 1;
                            w(counter) = stf(i).ray(j).weight(k);
                        end
                    end
                end
            end            
            
            %Set direct dose calculation and compute "dij"
            this.calcDoseDirect = true;
            dij = this.calcDose(ct,cst,stf);

            % calculate cubes; use uniform weights here, weighting with actual fluence 
            % already performed in dij construction
            resultGUI    = matRad_calcCubes(ones(dij.numOfBeams,1),dij);
            resultGUI.w  = w; 
        end

        function dij = calcDoseInfluence(this,ct,cst,stf)
            this.calcDoseDirect = false;
            dij = this.calcDose(ct,cst,stf);
        end
        function setDefaults(this)
            % future code for property validation on creation here
            matRad_cfg = MatRad_Config.instance();
            
            %Assign default parameters from MatRad_Config
            this.doseGrid.resolution    = matRad_cfg.propDoseCalc.defaultResolution;
            this.multScen = 'nomScen';
            this.selectVoxelsInScenarios = matRad_cfg.propDoseCalc.defaultSelectVoxelsInScenarios;
        end
    end
    
    methods(Access  =  protected)
        
        % method for setting and preparing the inition parameters for the 
        % dose calculation.
        % Should be called at the beginning of calcDose method.
        % Can be expanded or changed by overwriting this method and calling
        % the superclass method inside of it
        [dij,ct,cst,stf] = initDoseCalc(this,ct,cst,stf)   
        
        % method for finalizing the dose calculation (e.g. postprocessing
        % on dij or files
        function dij = finalizeDose(this,dij)

            matRad_cfg = MatRad_Config.instance();
            %Close Waitbar
            if any(ishandle(this.hWaitbar))
                delete(this.hWaitbar);
            end

            this.timers.full = toc(this.timers.full);
            
            matRad_cfg.dispInfo('Dose calculation finished in %g seconds!\n',this.timers.full);
        end

    
        function progressUpdate(this,pos,total)
            if nargin < 3
                pos = pos*1000;
                total=1000;
            end
            
            if pos ~= total && toc(this.lastProgressUpdate) < 1e-1
                return;
            end

            matRad_progress(pos,total);
            if any(ishandle(this.hWaitbar))
                waitbar(pos/total,this.hWaitbar);
            end
            
            this.lastProgressUpdate = tic;
        end
    end
    
    % Should be abstract methods but in order to satisfy the compatibility
    % with OCTAVE we can't use abstract methods. If OCTAVE at some point 
    % in the far future implements this feature this should be abstract again.
    methods (Access = protected) %(Abstract)                
        % the actual calculation method wich returns the final dij struct.
        % Needs to be implemented in non abstract subclasses. 
        % (Internal logic is often split into multiple methods in order to make the whole calculation more modular)
        function dij = calcDose(this,ct,cst,stf)
            throw(MException('MATLAB:class:AbstractMember','Abstract function calcDose of your DoseEngine needs to be implemented!'));
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
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('This is an Abstract Base class! Function needs to be called for instantiable subclasses!');
        end
        
        % static factory method to create/get correct dose engine from pln
        engine = getEngineFromPln(pln);

        % Machine Loader
        % Currently just uses the matRad function that asks for pln
        function machine = loadMachine(radiationMode,machineName)
            machine = matRad_loadMachine(struct('radiationMode',radiationMode,'machine',machineName));
        end
    end


end
