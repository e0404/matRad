classdef matRad_StfGeneratorBase < handle
% matRad_StfGeneratorBase: Abstract Superclass for Steering information 
%   generators. Steering information is used to guide the dose calculation
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

    properties (Constant, Abstract)
        name;                       %Descriptive Name
        shortName;                  %Short name for referencing
        possibleRadiationModes;     %Possible radiation modes for the respective StfGenerator
    end

    properties
        visMode = 0;                %Visualization Options
        addMargin = true;           %Add margins to target (for numerical stability and robustness)
        multScen;                   %Scenario Model
        bioParam;                   %Biological Model
        radiationMode;              %Radiation Mode
        machine;                    %Machine
    end

    properties (Access = protected)
        pln                         %matRad Plan struct
        cubeDim                     %Underlying CT dimension
        voxTargetWorldCoords        %Target voxels in world coordinates
        voiTarget                   %merged Target VOI cube
        ct                          %ct reference during generation
        cst                         %cst reference during generation
    end

    properties (Constant)
        isStfGenerator = true;    % const boolean for checking inheritance
    end

    methods 
        function this = matRad_StfGeneratorBase(pln)
            % Constructs standalone StfGenerator with or without pln

            this.setDefaults();
            if nargin == 1 && ~isempty(pln)
                this.assignPropertiesFromPln(pln);
            end
            
        end

        function setDefaults(this)
            % set default values from MatRad_Config

            matRad_cfg = MatRad_Config.instance();
            defaultPropStf = matRad_cfg.defaults.propStf;
            fields = fieldnames(defaultPropStf);
            for i = 1:numel(fields)
                fName = fields{i};
                if matRad_ispropCompat(this,fName)
                    try
                        this.(fName) = defaultPropStf.(fName);
                    catch
                        matRad_cfg.dispWarning('Could not assign default property %s',fName);
                    end
                end
            end
            this.machine = 'Generic';
        end

        function set.radiationMode(this,mode)
            % radiationMode setter

            if ~any(strcmp(mode,this.possibleRadiationModes))
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Radiation mode %s not supported by stf generator ''%s'' (%s)!',mode,this.name,class(this));
            end
            this.radiationMode = mode;
        end
        
        function warnDeprecatedProperty(this,oldProp,msg,newProp)
            matRad_cfg = MatRad_Config.instance();
            if nargin < 3 || isempty(msg)
                msg = '';
            end

            if nargin < 4
                dep2 = '';
            else
                dep2 = sprintf('Use Property ''%s'' instead!',newProp);
            end

            matRad_cfg.dispDeprecationWarning('Property ''%s'' of stf Generator ''%s'' is deprecated! %s%s',oldProp,this.name,msg,dep2);
        end

        function assignPropertiesFromPln(this,pln,warnWhenPropertyChanged)
            % Assign properties from pln.propStf to the stf generator

            matRad_cfg = MatRad_Config.instance();
            
            %Set Scenario Model
            if isfield(pln,'multScen')
                this.multScen = pln.multScen;
            end
            
            %Assign biological model
            if isfield(pln,'bioParam')
                this.bioParam = pln.bioParam;
            end

            %Take machine from pln
            if isfield(pln,'machine')
                this.machine = pln.machine;
            end
            
            if nargin < 3 || ~isscalar(warnWhenPropertyChanged) || ~islogical(warnWhenPropertyChanged)
                warnWhenPropertyChanged = false;
            end

            %Overwrite default properties within the stf generatorwith the 
            %ones given in the propStf struct
            if isfield(pln,'propStf') && isstruct(pln.propStf)
                plnStruct = pln.propStf; %get remaining fields
                if isfield(plnStruct,'type') && ~isempty(plnStruct.type) && ~any(strcmp(plnStruct.type,this.shortName))
                    matRad_cfg.dispWarning('Inconsistent stf generators given! pln asks for ''%s'', but you are using ''%s''!',plnStruct.type,this.shortName);
                end
                if isfield(plnStruct,'type')
                    plnStruct = rmfield(plnStruct, 'type'); % type field is no longer needed and would throw an exception
                end
            else
                plnStruct = struct();
            end

            fields = fieldnames(plnStruct);
            
            %Set up warning message
            if warnWhenPropertyChanged
                warningMsg = 'Property in stf generator overwritten from pln.propStf';
            else
                warningMsg = '';
            end

            % iterate over all fieldnames and try to set the
            % corresponding properties inside the stf generator
            if matRad_cfg.isOctave
                c2sWarningState = warning('off','Octave:classdef-to-struct');                
            end
            
            for i = 1:length(fields)
                try
                    field = fields{i};
                    if matRad_ispropCompat(this,field)
                        this.(field) = matRad_recursiveFieldAssignment(this.(field),plnStruct.(field),true,warningMsg);
                    else
                        matRad_cfg.dispWarning('Not able to assign property ''%s'' from pln.propStf to stf generator!',field);
                    end
                catch ME
                % catch exceptions when the stf generator has no 
                % properties which are defined in the struct.
                % When defining an engine with custom setter and getter
                % methods, custom exceptions can be caught here. Be
                % careful with Octave exceptions!
                    if ~isempty(warningMsg)
                        matRad_cfg = MatRad_Config.instance();
                        switch ME.identifier
                            case 'MATLAB:noPublicFieldForClass'
                                matRad_cfg.dispWarning('Not able to assign property from pln.propStf to stf generator: %s',ME.message);
                            otherwise
                                matRad_cfg.dispWarning('Problem while setting up stf generator from struct:%s %s',field,ME.message);
                        end
                    end
                end
            end
            
            if matRad_cfg.isOctave
                warning(c2sWarningState.state,'Octave:classdef-to-struct');                
            end
        end
    end

    methods 

        function stf = generate(this, ct, cst)  
            % Generate steering information for the given ct and cst
            % This is a base class function performing the following tasks
            % 1. Checking the input
            % 2. Initializing the patient geometry (protected properties)
            % 3. Generating the source information (and thus the "stf")

            % Instance of MatRad_Config class
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispInfo('matRad: Generating stf struct with generator ''%s''... ',this.name);

            if isempty(this.multScen)
                this.multScen = matRad_NominalScenario(ct);
            end

            % get machine
            if ~isstruct(this.machine)
                try
                    this.machine = matRad_loadMachine(struct('radiationMode',this.radiationMode,'machine',this.machine));
                catch
                    matRad_cfg.dispError('Could not find the following machine file: %s_%s',this.radiationMode,this.machine);
                end
            end

            this.ct = ct;
            this.cst = cst;

            this.initializePatientGeometry();
            stf = this.generateSourceGeometry();
        end
    end

    methods (Access = protected)

        function initializePatientGeometry(this)
            % Basic Initialization of the Patient Geometry

            matRad_cfg = MatRad_Config.instance();
            
            % Initialize patient geometry
            V = [];
            %ct = matRad_calcWaterEqD(ct,this.radiationMode);  

            isTarget = cellfun(@(voiType) isequal(voiType, 'TARGET'), this.cst(:,3));
            if ~any(isTarget)
                matRad_cfg.dispError('No target found in cst. Please designate at least one VOI as ''TARGET''!');
            end

            hasObjective = ~cellfun(@isempty, this.cst(:,6));
            useTargetForBixelPlacement = isTarget & hasObjective;

            if ~any(useTargetForBixelPlacement)
                matRad_cfg.dispWarning('No Objectives / Constraints assigned to targets. All targets will be considered for Bixel placement!');
                useTargetForBixelPlacement(isTarget) = true;
            end

            % Now add all used target voxels to the voxel list
            for i = 1:size(this.cst, 1)
                if useTargetForBixelPlacement(i)
                    V = [V; this.cst{i,4}{1}];  
                end
            end

            % Remove double voxels
            V = unique(V);
            % generate voi cube for targets
            this.voiTarget = zeros(this.ct.cubeDim);
            this.voiTarget(V) = 1;

            % Margin info
            if this.addMargin
                pbMargin = this.getPbMargin();

                % Assumption for range uncertainty
                assumeRangeMargin = this.multScen.maxAbsRangeShift + this.multScen.maxRelRangeShift + pbMargin;

                % add margin - account for voxel resolution, the maximum shift scenario and the current bixel width.
                margin.x = max([this.ct.resolution.x max(abs(this.multScen.isoShift(:,1)) + assumeRangeMargin)]);
                margin.y = max([this.ct.resolution.y max(abs(this.multScen.isoShift(:,2)) + assumeRangeMargin)]);
                margin.z = max([this.ct.resolution.z max(abs(this.multScen.isoShift(:,3)) + assumeRangeMargin)]);

                this.voiTarget = matRad_addMargin(this.voiTarget, this.cst, this.ct.resolution, margin, true);
                V = find(this.voiTarget > 0);
            end

            % throw error message if no target is found
            if isempty(V)
                matRad_cfg.dispError('Could not find target.');
            end

            % get world coordinate system
            this.ct = matRad_getWorldAxes(this.ct);

            % Convert linear indices to 3D voxel coordinates
            this.voxTargetWorldCoords = matRad_cubeIndex2worldCoords(V, this.ct);

            % take only voxels inside patient
            V = [this.cst{:,4}];
            V = unique(vertcat(V{:}));

            % ignore densities outside of contours
            eraseCtDensMask = ones(prod(this.ct.cubeDim), 1);
            eraseCtDensMask(V) = 0;
            for i = 1:this.ct.numOfCtScen
                this.ct.cube{i}(eraseCtDensMask == 1) = 0;
            end

            
        end

        function pbMargin = getPbMargin(this)
            % Get the geometrical margin to add to the target (e.g. for safe spot placement or because of robustness)
            pbMargin = 0;
        end
    end

    methods (Access = protected)
        function stf = generateSourceGeometry(this)
            % the actual calculation method which returns the final stf struct.
            % Needs to be implemented in non-abstract subclasses.
            % (Internal logic is often split into multiple methods in order to
            % make the whole calculation more modular)
            throw(MException('MATLAB:class:AbstractMember','Abstract function generateSourceGeometry of your StfGenerator needs to be implemented!'));
        end
    end
end
