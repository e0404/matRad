classdef (Abstract) matRad_StfGeneratorBase < handle
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

    properties (Constant)
        isStfGenerator = true;      % const boolean for inheritance quick check
    end

    properties (Constant, Abstract)
        name;                       %Descriptive Name
        shortName;                  %Short name for referencing
        possibleRadiationModes;     %Possible radiation modes for the respective StfGenerator
    end

    properties (Access = public)
        visMode = 0;                %Visualization Options
        addMargin = true;           %Add margins to target (for numerical stability and robustness)
        multScen;                   %Scenario Model
        bioModel;                   %Biological Model
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
            
            %Must haves in pln struct
            %Set/validate radiation Mode
            if ~isfield(pln,'radiationMode') && isempty(this.radiationMode)
                matRad_cfg.dispError('No radiation mode specified in pln struct!');
            else
                this.radiationMode = pln.radiationMode;
            end

            %Take machine from pln
            if ~isfield(pln,'machine') && isempty(this.machine)
                matRad_cfg.dispError('No machine specified in pln struct!');
            else
                this.machine = pln.machine;
            end
            
            % Defaults if not provided
            %Set Scenario Model
            if isfield(pln,'multScen')
                this.multScen = pln.multScen;
            end

            %Assign biological model
            if isfield(pln,'bioModel')
                try
                    this.bioModel = matRad_BiologicalModel.validate(pln.bioModel, pln.radiationMode);
                catch ME
                    %Steering generation usually works independent of biological model, so we warn and set a dummy model
                    matRad_cfg.dispWarning('Biological model inconsistent with chosen machine / radiation mode, biological dose calculation will probably not work: %s',ME.message);
                    this.bioModel = matRad_EmptyBiologicalModel();
                end
            else
                this.bioModel = matRad_EmptyBiologicalModel();
            end

            if nargin < 3 || ~isscalar(warnWhenPropertyChanged) || ~islogical(warnWhenPropertyChanged)
                warnWhenPropertyChanged = false;
            end

            %Overwrite default properties within the stf generatorwith the
            %ones given in the propStf struct
            if isfield(pln,'propStf') && isstruct(pln.propStf)
                plnStruct = pln.propStf; %get remaining fields
                if isfield(plnStruct,'generator') && ~isempty(plnStruct.generator) && ~any(strcmp(plnStruct.generator,this.shortName))
                    matRad_cfg.dispWarning('Inconsistent stf generators given! pln asks for ''%s'', but you are using ''%s''!',plnStruct.generator,this.shortName);
                end
                if isfield(plnStruct,'generator')
                    plnStruct = rmfield(plnStruct, 'generator'); % generator field is no longer needed and would throw an exception
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

    methods %(Abstract)

        function stf = generate(this, ct, cst)
            % Generate steering information for the given ct and cst
            % This is a base class function performing the following tasks
            % 1. Checking the input
            % 2. Initializing the patient geometry (protected properties)
            % 3. Generating the source information (and thus the "stf")

            % Instance of MatRad_Config class
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispInfo('matRad: Generating stf struct with generator ''%s''... ',this.name);

            this.ct = ct;
            this.cst = cst;

            this.initialize();
            this.createPatientGeometry();
            stf = this.generateSourceGeometry();
        end
    end

    methods (Access = protected)
        function initialize(this)
            %Do nothing
            matRad_cfg = MatRad_Config.instance();
            % get machine
            if ~isstruct(this.machine)
                try
                    this.machine = matRad_loadMachine(struct('radiationMode',this.radiationMode,'machine',this.machine));
                catch
                    matRad_cfg.dispError('Could not find the following machine file: %s_%s',this.radiationMode,this.machine);
                end
            end

            %Make sure the coordinate system is initialized
            this.ct = matRad_getWorldAxes(this.ct);

            %Make sure we have a  valid scenario model
            if isempty(this.multScen)
                this.multScen = matRad_NominalScenario(this.ct);
            end
            
            %create / validate scenario model
            this.multScen = matRad_ScenarioModel.create(this.multScen,this.ct);
        end

        function createPatientGeometry(this)
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

    methods (Static)
        function generator = getGeneratorFromPln(pln, warnDefault)
            %GETENGINE Summary of this function goes here
            %   Detailed explanation goes here

            if nargin < 2
                warnDefault = true;
            end

            matRad_cfg = MatRad_Config.instance();

            generator = [];

            initDefaultGenerator = false;
            %get all available engines for given pln struct, could be done conditional
            classList = matRad_StfGeneratorBase.getAvailableGenerators(pln);

            % Check for a valid engine, and if the given engine isn't valid set boolean
            % to initiliaze default engine at the end of this function
            if isfield(pln,'propStf') && isa(pln.propStf, mfilename('class'))
                generator = pln.propStf;
            elseif isfield(pln,'propStf') && isstruct(pln.propStf) && isfield(pln.propStf,'generator')
                if ischar(pln.propStf.generator) || isstring(pln.propStf.generator)
                    matchGenerators = strcmpi({classList(:).shortName},pln.propStf.generator);
                    if any(matchGenerators)
                        %instantiate engine
                        generatorHandle = classList(matchGenerators).handle;
                        generator = generatorHandle(pln);
                    else
                        initDefaultGenerator = true;
                    end
                else
                    initDefaultGenerator = true;
                    matRad_cfg.dispWarning('pln.propStf.generator field not valid!');
                end
            else
                initDefaultGenerator = true;
            end

            % trying to use a default engine which fits
            % the given radiation mode, when no valid engine was defined.
            % Default Engines are defined in matRad_Config.
            if initDefaultGenerator
                matchGenerators = ismember({classList(:).shortName},matRad_cfg.defaults.propStf.generator);
                if any(matchGenerators)
                    generatorHandle = classList(matchGenerators).handle;

                    % unlikely event that multiple engines fit just take the first
                    if length(generatorHandle) > 1
                        generatorHandle = generatorHandle{1};
                    end
                    generator = generatorHandle(pln);
                    if warnDefault
                        matRad_cfg.dispWarning('Using default stf generator %s!', generator.name);
                    end
                elseif ~isempty(classList)
                    generatorHandle = classList(1).handle;
                    generator = generatorHandle(pln);
                    matRad_cfg.dispWarning('Default stf generator not available! Using %s.', generator.name);
                else
                    matRad_cfg.dispError('Default stf generator not found!');
                end
            end

            if isempty(generator)
                matRad_cfg.dispError('No suitable stf generator found!');
            end

        end

        function classList = getAvailableGenerators(pln,optionalPaths)
            % Returns a list of names and coresponding handle for stf
            % generators. Returns all stf generators when no arg is
            %   given. If no generators are found return gonna be empty.
            %
            % call:
            %   classList = matRad_StfGeneratorBase.getAvailableGenerators(pln,optional_path)
            %
            % input:
            %   pln:            containing proposed generator and machine file informations
            %   optionalPath:   cell array of other folders to search in
            %
            % returns:
            %   classList: struct array containing name, shortName, className, and
            %   handle for construction of the Generator

            matRad_cfg = MatRad_Config.instance();

            %Parse inputs
            if nargin < 2
                optionalPaths = {fileparts(mfilename("fullpath"))};
            else
                if ~(iscellstr(optionalPaths) && all(optionalPaths))
                    matRad_cfg.dispError('Invalid path array!');
                end

                optionalPaths = horzcat(fileparts(mfilename("fullpath")),optionalPaths);
            end

            if nargin < 1
                pln = [];
            else
                if ~(isstruct(pln) || isempty(pln))
                    matRad_cfg.dispError('Invalid pln!');
                end
            end

            %Get available, valid classes through call to matRad helper function
            %for finding subclasses
            persistent allAvailableStfGenerators lastOptionalPaths
            if isempty(allAvailableStfGenerators) || (~isempty(lastOptionalPaths) && ~isequal(lastOptionalPaths, optionalPaths))
                lastOptionalPaths = optionalPaths;
                allAvailableStfGenerators = matRad_findSubclasses(mfilename('class'),'folders',optionalPaths,'includeAbstract',false);
            end

            availableStfGenerators = allAvailableStfGenerators;

            %Now filter for pln
            ix = [];

            if nargin >= 1 && ~isempty(pln)
                machine = matRad_loadMachine(pln);
                machineMode = machine.meta.radiationMode;

                for cIx = 1:length(availableStfGenerators)
                    mc = availableStfGenerators{cIx};
                    availabilityFuncStr = [mc.Name '.isAvailable'];
                    %availabilityFunc = str2func(availabilityFuncStr); %str2func  does not seem to work on static class functions in Octave 5.2.0
                    try
                        %available = availabilityFunc(pln,machine);
                        available = eval([availabilityFuncStr '(pln,machine)']);
                    catch
                        available = false;
                        mpList = mc.PropertyList;
                        if matRad_cfg.isMatlab
                            loc = find(arrayfun(@(x) strcmp('possibleRadiationModes',x.Name),mpList));
                            propValue = mpList(loc).DefaultValue;
                        else
                            loc = find(cellfun(@(x) strcmp('possibleRadiationModes',x.Name),mpList));
                            propValue = mpList{loc}.DefaultValue;
                        end

                        if any(strcmp(propValue, pln.radiationMode))
                            % get radiation mode from the in pln proposed basedata machine file
                            % add current class to return lists if the
                            % radiation mode is compatible
                            if(any(strcmp(propValue, machineMode)))
                                available = true;

                            end
                        end
                    end
                    if available
                        ix = [ix cIx];
                    end
                end

                availableStfGenerators = availableStfGenerators(ix);
            end

            classList = matRad_identifyClassesByConstantProperties(availableStfGenerators,'shortName','defaults',matRad_cfg.defaults.propStf.generator,'additionalPropertyNames',{'name'});

        end

        function [available,msg] = isAvailable(pln,machine)
            % return a boolean if the generator is is available for the given pln
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

        % Machine Loader
        % Currently just uses the matRad function that asks for pln
        function machine = loadMachine(radiationMode,machineName)
            machine = matRad_loadMachine(struct('radiationMode',radiationMode,'machine',machineName));
        end
    end
end
