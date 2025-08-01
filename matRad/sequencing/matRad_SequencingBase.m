classdef matRad_SequencingBase < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Constant)
        isSequencer = true;      % const boolean for inheritance quick check
    end

    properties (Constant, Abstract)
        name;                       %Descriptive Name
        shortName;                  %Short name for referencing
        possibleRadiationModes;     %Possible radiation modes for the respective Sequencer
    end

    properties (Access = public)
        radiationMode;              %Radiation Mode
        visMode = 0;                % vis bool
    end

    methods
        function this = matRad_SequencingBase(pln)
            % Constructs standalone sequencer with or without pln

            this.setDefaults();
            if nargin == 1 && ~isempty(pln)
                this.assignPropertiesFromPln(pln);
            end

        end

        function setDefaults(this)
            % set default values from MatRad_Config

            matRad_cfg = MatRad_Config.instance();
            defaultPropSeq = matRad_cfg.defaults.propSeq;
            fields = fieldnames(defaultPropSeq);
            for i = 1:numel(fields)
                fName = fields{i};
                if matRad_ispropCompat(this,fName)
                    try
                        this.(fName) = defaultPropSeq.(fName);
                    catch
                        matRad_cfg.dispWarning('Could not assign default property %s',fName);
                    end
                end
            end
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

            matRad_cfg.dispDeprecationWarning('Property ''%s'' of sequencer ''%s'' is deprecated! %s%s',oldProp,this.name,msg,dep2);
        end
        
        function assignPropertiesFromPln(this,pln,warnWhenPropertyChanged)
            % Assign properties from pln.propSeq to the sequencer

            matRad_cfg = MatRad_Config.instance();
            
            %Must haves in pln struct
            %Set/validate radiation Mode
            if ~isfield(pln,'radiationMode') && isempty(this.radiationMode)
                matRad_cfg.dispError('No radiation mode specified in pln struct!');
            else
                this.radiationMode = pln.radiationMode;
            end

            if nargin < 3 || ~isscalar(warnWhenPropertyChanged) || ~islogical(warnWhenPropertyChanged)
                warnWhenPropertyChanged = false;
            end

            %Overwrite default properties within the sequencer with the
            %ones given in the propSeq struct
            if isfield(pln,'propSeq') && isstruct(pln.propSeq)
                plnStruct = pln.propSeq; %get remaining fields
                if isfield(plnStruct,'sequencer') && ~isempty(plnStruct.sequencer) && ~any(strcmp(plnStruct.sequencer,this.shortName))
                    matRad_cfg.dispWarning('Inconsistent sequencers given! pln asks for ''%s'', but you are using ''%s''!',plnStruct.sequencer,this.shortName);
                end
                if isfield(plnStruct,'sequencer')
                    plnStruct = rmfield(plnStruct, 'sequencer'); % sequencer field is no longer needed and would throw an exception
                end
            else
                plnStruct = struct();
            end

            fields = fieldnames(plnStruct);

            %Set up warning message
            if warnWhenPropertyChanged
                warningMsg = 'Property in sequencer overwritten from pln.propSeq';
            else
                warningMsg = '';
            end

            % iterate over all fieldnames and try to set the
            % corresponding properties inside the sequencer
            if matRad_cfg.isOctave
                c2sWarningState = warning('off','Octave:classdef-to-struct');
            end

            for i = 1:length(fields)
                try
                    field = fields{i};
                    if matRad_ispropCompat(this,field)
                        this.(field) = matRad_recursiveFieldAssignment(this.(field),plnStruct.(field),true,warningMsg);
                    else
                        matRad_cfg.dispWarning('Not able to assign property ''%s'' from pln.propSeq to sequencer',field);
                    end
                catch ME
                    % catch exceptions when the sequencing has no
                    % properties which are defined in the struct.
                    % When defining an engine with custom setter and getter
                    % methods, custom exceptions can be caught here. Be
                    % careful with Octave exceptions!
                    if ~isempty(warningMsg)
                        matRad_cfg = MatRad_Config.instance();
                        switch ME.identifier
                            case 'MATLAB:noPublicFieldForClass'
                                matRad_cfg.dispWarning('Not able to assign property from pln.propSeq to sequencing: %s',ME.message);
                            otherwise
                                matRad_cfg.dispWarning('Problem while setting up sequencing from struct:%s %s',field,ME.message);
                        end
                    end
                end
            end

            if matRad_cfg.isOctave
                warning(c2sWarningState.state,'Octave:classdef-to-struct');
            end
        end

    end
    methods (Static)
        function sequencer = getSequencerFromPln(pln, warnDefault)
            %GETENGINE Summary of this function goes here
            %   Detailed explanation goes here

            if nargin < 2
                warnDefault = true;
            end

            matRad_cfg = MatRad_Config.instance();

            sequencer = [];

            initDefaultSequencer = false;
            %get all available Sequencers for given pln struct, could be done conditional
            classList = matRad_SequencingBase.getAvailableSequencers(pln);

            % Check for a valid engine, and if the given engine isn't valid set boolean
            % to initiliaze default engine at the end of this function
            if isfield(pln,'propSeq') && isa(pln.propSeq, mfilename('class'))
                sequencer = pln.propSeq;
            elseif isfield(pln,'propSeq') && isstruct(pln.propSeq) && isfield(pln.propSeq,'sequencer')
                if ischar(pln.propSeq.sequencer) || isstring(pln.propSeq.sequencer)
                    matchSequencers = strcmpi({classList(:).shortName},pln.propSeq.sequencer);
                    if any(matchSequencers)
                        %instantiate engine
                        sequencerHandle = classList(matchSequencers).handle;
                        sequencer = sequencerHandle(pln);
                    else
                        initDefaultSequencer = true;
                    end
                else
                    initDefaultSequencer = true;
                    matRad_cfg.dispWarning('pln.propSeq.sequencer field not valid!');
                end
            else
                initDefaultSequencer = true;
            end

            % trying to use a default engine which fits
            % the given radiation mode, when no valid engine was defined.
            % Default Engines are defined in matRad_Config.
            if initDefaultSequencer
                matchSequencers = ismember({classList(:).shortName},matRad_cfg.defaults.propSeq.sequencer);
                if any(matchSequencers)
                    sequencerHandle = classList(matchSequencers).handle;

                    % unlikely event that multiple engines fit just take the first
                    if length(sequencerHandle) > 1
                        sequencerHandle = sequencerHandle{1};
                    end
                    sequencer = sequencerHandle(pln);
                    if warnDefault
                        matRad_cfg.dispWarning('Using default sequencer %s!', sequencer.name);
                    end
                elseif ~isempty(classList)
                    sequencerHandle = classList(1).handle;
                    sequencer = sequencerHandle(pln);
                    matRad_cfg.dispWarning('Default sequencer not available! Using %s.', sequencer.name);
                else
                    matRad_cfg.dispError('Default sequencer not found!');
                end
            end

            if isempty(sequencer)
                matRad_cfg.dispError('No suitable sequencer found!');
            end

        end

        function classList = getAvailableSequencers(pln,optionalPaths)
     
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
            persistent allAvailableSequencers lastOptionalPaths
            
            %First we do a sanity check if persistently stored metaclasses are valid
            if ~matRad_cfg.isOctave && ~isempty(allAvailableSequencers) && ~all(cellfun(@isvalid,allAvailableSequencers))
                matRad_cfg.dispWarning('Found invalid Sequencing Sequencers, updating cache.');
                allAvailableSequencers = [];
            end

            if isempty(allAvailableSequencers) || (~isempty(lastOptionalPaths) && ~isequal(lastOptionalPaths, optionalPaths))
                lastOptionalPaths = optionalPaths;
                allAvailableSequencers = matRad_findSubclasses(mfilename('class'),'folders',optionalPaths,'includeAbstract',false);
            end

           availableSequencers = allAvailableSequencers;

            %Now filter for pln
            ix = [];

            if nargin >= 1 && ~isempty(pln)
                machine = matRad_loadMachine(pln);
                machineMode = machine.meta.radiationMode;

                for cIx = 1:length(availableSequencers)
                    mc = availableSequencers{cIx};
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

                availableSequencers = availableSequencers (ix);
            end

            classList = matRad_identifyClassesByConstantProperties(availableSequencers,'shortName','defaults',matRad_cfg.defaults.propSeq.sequencer,'additionalPropertyNames',{'name'});

        end

        function [available,msg] = isAvailable(pln,machine)

            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('This is an Abstract Base class! Function needs to be called for instantiable subclasses!');
        end
    end

end