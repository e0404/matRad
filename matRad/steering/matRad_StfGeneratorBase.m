classdef matRad_StfGeneratorBase < handle

    properties (Constant, Abstract)
        name;
        shortName;
        possibleRadiationModes;
    end

    properties
        visMode = 0;
        addMargin = true;
        multScen;
        bioParam;
        radiationMode;
        machine;
    end

    properties (Access = protected)
        pln
        cubeDim
        voxTargetWorldCoords
        voiTarget
    end

    properties (Constant)
        isStfGenerator = true;    % const boolean for checking inheritance
    end

    methods 
        function this = matRad_StfGeneratorBase(pln)
            this.setDefaults();
            if nargin == 1 && ~isempty(pln)
                this.assignPropertiesFromPln(pln);
            end
            this.radiationMode = 'brachy';
        end

        function setDefaults(this)
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

            this.initializePatientGeometry(ct, cst);
            stf = this.generateSourceGeometry(ct, cst);
        end
    end

    methods (Access = protected)

        function initializePatientGeometry(this, ct, cst)
            matRad_cfg = MatRad_Config.instance();
            
            % Initialize patient geometry
            V = [];
            %ct = matRad_calcWaterEqD(ct,this.radiationMode);  

            isTarget = cellfun(@(voiType) isequal(voiType, 'TARGET'), cst(:,3));
            if ~any(isTarget)
                matRad_cfg.dispError('No target found in cst. Please designate at least one VOI as ''TARGET''!');
            end

            hasObjective = ~cellfun(@isempty, cst(:,6));
            useTargetForBixelPlacement = isTarget & hasObjective;

            if ~any(useTargetForBixelPlacement)
                matRad_cfg.dispWarning('No Objectives / Constraints assigned to targets. All targets will be considered for Bixel placement!');
                useTargetForBixelPlacement(isTarget) = true;
            end

            % Now add all used target voxels to the voxel list
            for i = 1:size(cst, 1)
                if useTargetForBixelPlacement(i)
                    V = [V; cst{i,4}{1}];  
                end
            end

            % Remove double voxels
            V = unique(V);
            % generate voi cube for targets
            this.voiTarget = zeros(ct.cubeDim);
            this.voiTarget(V) = 1;

            % Margin info
            if this.addMargin
                pbMargin = this.getPbMargin();

                % Assumption for range uncertainty
                assumeRangeMargin = this.multScen.maxAbsRangeShift + this.multScen.maxRelRangeShift + pbMargin;

                % add margin - account for voxel resolution, the maximum shift scenario and the current bixel width.
                margin.x = max([ct.resolution.x max(abs(this.multScen.isoShift(:,1)) + assumeRangeMargin)]);
                margin.y = max([ct.resolution.y max(abs(this.multScen.isoShift(:,2)) + assumeRangeMargin)]);
                margin.z = max([ct.resolution.z max(abs(this.multScen.isoShift(:,3)) + assumeRangeMargin)]);

                this.voiTarget = matRad_addMargin(this.voiTarget, cst, ct.resolution, margin, true);
                V = find(this.voiTarget > 0);
            end

            % throw error message if no target is found
            if isempty(V)
                matRad_cfg.dispError('Could not find target.');
            end

            % get world coordinate system
            ct = matRad_getWorldAxes(ct);

            % Convert linear indices to 3D voxel coordinates
            this.voxTargetWorldCoords = matRad_cubeIndex2worldCoords(V, ct);

            % take only voxels inside patient
            V = [cst{:,4}];
            V = unique(vertcat(V{:}));

            % ignore densities outside of contours
            eraseCtDensMask = ones(prod(ct.cubeDim), 1);
            eraseCtDensMask(V) = 0;
            for i = 1:ct.numOfCtScen
                ct.cube{i}(eraseCtDensMask == 1) = 0;
            end

            
        end

        function pbMargin = getPbMargin(this)
            pbMargin = 0;
        end
    end

    methods (Access = protected)
        % the actual calculation method which returns the final stf struct.
        % Needs to be implemented in non-abstract subclasses.
        % (Internal logic is often split into multiple methods in order to
        % make the whole calculation more modular)
        function stf = generateSourceGeometry(this, ct, cst)
            throw(MException('MATLAB:class:AbstractMember','Abstract function generateSourceGeometry of your StfGenerator needs to be implemented!'));
        end
    end
end
