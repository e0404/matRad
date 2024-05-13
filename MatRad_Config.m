classdef MatRad_Config < handle
    % MatRad_Config MatRad Configuration class
    % This class is used globally through Matlab to handle default values and
    % logging and is declared as global matRad_cfg.
    % Usage:
    %    matRad_cfg = MatRad_Config.instance();
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2019 the matRad development team.
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

        %Logging
        logLevel = 4; %1 = only Errors, 2 = with Warnings, 3 = Info output, 4 = deprecation warnings, 5 = debug information
        keepLog = false; %Stores the full log in memory
        writeLog = false; %Writes the log to a file on-the-fly

        %Default Properties
        propDoseCalc;
        propOpt;
        propMC;
        propStf;

        defaults;

        %Disable GUI
        disableGUI = false;

        devMode = false;
        eduMode = false;
        
        gui;
    end

    properties (SetAccess = private)
        messageLog = {};
        logFileHandle;

        %For storing the Environment & its version
        env;
        envVersion;
        isOctave; %Helper bool to check for Octave
        isMatlab; %Helper bool to check for Matlab
        matRad_version; %MatRad version string

        matRadRoot; %Path to matRadRoot
    end

    methods (Access = private)
        function obj = MatRad_Config()
            %MatRad_Config Constructs an instance of this class.
            %  The configuration is implemented as a singleton and used globally
            %  Therefore its constructor is private
            %  For instantiation, use the static MatRad_Config.instance();
            
            %Set Path
            obj.matRadRoot = fileparts(mfilename('fullpath'));
            addpath(genpath(obj.matRadRoot));

            %Set Version
            obj.getEnvironment();
            obj.matRad_version = matRad_version();

            %Configure Environment
            obj.configureEnvironment();

            %Just to catch people messing with the properties in the file
            if ~isempty(obj.writeLog) && obj.writeLog
                logFile = [obj.matRadRoot filesep 'matRad.log'];
                obj.logFileHandle = fopen(logFile,'a');
            end

            %Call the reset function for remaining inatialization
            obj.reset();
        end

        function delete(~)
            %might not be desired by users
            %rmpath(genpath(matRad_cfg.matRadRoot));
        end

        function displayToConsole(obj,type,formatSpec,varargin)
            %displayToConsole lowest-level logging function for matRad.
            %   Display to console will be called from the public wrapper
            %   functions dispError, dispWarning, dispInfo, dispDebug
            %
            %  input
            %    type:			type of the log information.
            %                   Needs to be one of 'error', 'warning', 'info' or 'debug'.
            %    formatSpec: 	string to print using format specifications similar to fprintf
            %    varargin:   	variables according to formatSpec

            if nargin < 4
                forwardArgs = {formatSpec};
            else
                forwardArgs = [{formatSpec},varargin(:)'];
            end

            if obj.keepLog
                obj.messageLog{end+1,1} = upper(type);
                obj.messageLog{end,2} = sprintf(forwardArgs{:});
            end

            switch type
                case{'info'}
                    if obj.logLevel >= 3
                        fprintf(forwardArgs{:});
                    end
                case{'debug'}
                    if obj.logLevel >= 5
                        forwardArgs{1} = ['DEBUG: ' forwardArgs{1}];
                        fprintf(forwardArgs{:});
                    end
                case{'dep'}
                    if obj.logLevel >= 4
                        forwardArgs{1} = ['DEPRECATION WARNING: ' forwardArgs{1}];
                        warning(forwardArgs{:});
                    end
                case{'warning'}
                    if obj.logLevel >= 2
                        warning(forwardArgs{:});
                    end
                case {'error'}
                    if obj.logLevel >= 1
                        %We create an error structure to later clean the
                        %stack trace from the last two files/lines (i.e.,
                        %this function / file)

                        err.message = sprintf(forwardArgs{:});
                        err.identifier = 'matRad:Error';
                        %err.stack = dbstack(2);
                        error(err);

                    end
                otherwise
                    error('Log type %s not defined!',type);
            end

            if obj.writeLog
                fprintf(obj.logFileHandle,forwardArgs{:});
            end
        end
    end

    methods
        function reset(obj)
            %Set all default properties for matRad's computations
            obj.setDefaultProperties();
            obj.setDefaultGUIProperties();
        end

        function setDefaultProperties(obj)
            %setDefaultProperties set matRad's default computation
            %   properties
            %  input
            
            %Geometric stf options
            obj.propStf.defaultLongitudinalSpotSpacing = 2;
            obj.propStf.defaultAddMargin = true; %expand target for beamlet finding
            obj.propStf.defaultBixelWidth = 5;
          
            %Dose Calculation Options
            obj.propDoseCalc.defaultDoseEngines = {'SVDPB','HongPB'}; %Names for default engines used when no other is given
            obj.propDoseCalc.defaultResolution = struct('x',3,'y',3,'z',3); %[mm]
            obj.propDoseCalc.defaultDosimetricLateralCutOff = 0.995; %[rel.]
            obj.propDoseCalc.defaultGeometricLateralCutOff = 50; %[mm]
            obj.propDoseCalc.defaultKernelCutOff = Inf; %[mm]
            obj.propDoseCalc.defaultSsdDensityThreshold = 0.05; %[rel.]
            obj.propDoseCalc.defaultUseGivenEqDensityCube = false; %Use the given density cube ct.cube and omit conversion from cubeHU.
            obj.propDoseCalc.defaultIgnoreOutsideDensities = true; %Ignore densities outside of cst contours
            obj.propDoseCalc.defaultUseCustomPrimaryPhotonFluence = false; %Use a custom primary photon fluence
            obj.propDoseCalc.defaultCalcLET = true; %calculate LETs for particles
            obj.propDoseCalc.defaultSelectVoxelsInScenarios = 'all';
            obj.propDoseCalc.defaultAirOffsetCorrection = true;

            % default properties for fine sampling calculation
            obj.propDoseCalc.defaultFineSamplingProperties.sigmaSub = 1;
            obj.propDoseCalc.defaultFineSamplingProperties.N = 2;
            obj.propDoseCalc.defaultFineSamplingProperties.method = 'fitCircle';

            %Monte Carlo options
            obj.propDoseCalc.defaultNumHistoriesPerBeamlet = 2e4;
            obj.propDoseCalc.defaultNumHistoriesDirect = 1e6;
            obj.propDoseCalc.defaultOutputMCvariance = true;

            %deprecated monte carlo options
            obj.propMC.ompMC_defaultHistories = 1e6;
            obj.propMC.ompMC_defaultOutputVariance = false;

            % Set default histories for MonteCarlo here if necessary
            %             obj.propMC.defaultNumHistories = 100;

            obj.propMC.default_photon_engine = 'matRad_OmpConfig';
            %             obj.propMC.default_photon_engine = 'matRad_TopasConfig';
            obj.propMC.default_proton_engine = 'matRad_MCsquareConfig';
            obj.propMC.default_carbon_engine = 'matRad_TopasConfig';

            % Default settings for TOPAS
            obj.propMC.default_beamProfile_particles = 'biGaussian';
            obj.propMC.default_beamProfile_photons = 'uniform';
            obj.propMC.defaultExternalCalculation = false;
            obj.propMC.defaultCalcDij = false;


            %Default Optimization Options
            obj.propOpt.defaultMaxIter = 10000;
            obj.propOpt.defaultRunDAO = 0;
            obj.propOpt.defaultRunSequencing = 0;
            obj.propOpt.defaultClearUnusedVoxels = false;
            


            obj.disableGUI = false;
            
            obj.defaults.samplingScenarios = 25;

            obj.devMode = false;
            obj.eduMode = false;

        end

        %%For testing
        function setDefaultPropertiesForTesting(obj)
            %setDefaultPropertiesForTesting sets matRad's default
            %properties during testing to reduce computational load

            obj.setDefaultProperties();

            obj.logLevel = 1; %Omit output except errors
            
            %Geomteric stf options
            obj.propStf.defaultLongitudinalSpotSpacing = 20;
            obj.propStf.defaultBixelWidth = 20;
            
            %Dose Calculation Options
            obj.propDoseCalc.defaultResolution = struct('x',5,'y',6,'z',7); %[mm]
            obj.propDoseCalc.defaultGeometricLateralCutOff = 20;
            obj.propDoseCalc.defaultDosimetricLateralCutOff = 0.8;
            obj.propDoseCalc.defaultKernelCutOff = 20; %[mm]
            obj.propDoseCalc.defaultSsdDensityThreshold = 0.05;
                      
            %Monte Carlo options
            obj.propDoseCalc.defaultNumHistoriesPerBeamlet = 100;
            obj.propDoseCalc.defaultNumHistoriesDirect = 100;
            obj.propDoseCalc.defaultOutputMCvariance = true;
            
            %Deprecated Monte Carlo options
            obj.propMC.ompMC_defaultHistories = 100;
            obj.propMC.ompMC_defaultOutputVariance = true;

            % Set default histories for MonteCarlo
            obj.propMC.defaultNumHistories = 100;

            % default optimization options
            obj.propOpt.defaultMaxIter = 10;

            obj.defaults.samplingScenarios = 2;

            obj.disableGUI = true;

            obj.devMode = true;
            obj.eduMode = false;
        end  
        
        %%for edu mode
        function setDefaultPropertiesForEduMode(obj)
            
            obj.setDefaultProperties();

            obj.logLevel = 1;
            
            obj.propStf.defaultLongitudinalSpotSpacing = 3;
            obj.propStf.defaultAddMargin = true; %expand target for beamlet finding
            obj.propStf.defaultBixelWidth = 5;
            
            obj.propDoseCalc.defaultResolution = struct('x',4,'y',4,'z',4); %[mm]
            obj.propDoseCalc.defaultLateralCutOff = 0.975; %[rel.]
            obj.propDoseCalc.defaultGeometricCutOff = 50; %[mm]
            obj.propDoseCalc.defaultSsdDensityThreshold = 0.05; %[rel.]
            obj.propDoseCalc.defaultUseGivenEqDensityCube = false; %Use the given density cube ct.cube and omit conversion from cubeHU.
            obj.propDoseCalc.defaultIgnoreOutsideDensities = true; %Ignore densities outside of cst contours
            obj.propDoseCalc.defaultUseCustomPrimaryPhotonFluence = false; %Use a custom primary photon fluence
            
            obj.propOpt.defaultMaxIter = 500;
            
            obj.propMC.ompMC_defaultHistories = 1e4;
            obj.propMC.ompMC_defaultOutputVariance = false;
            obj.propMC.MCsquare_defaultHistories = 1e4;
            obj.propMC.direct_defaultHistories = 1e4;
            
            obj.disableGUI = false;
            
            obj.devMode = false;
            obj.eduMode = true;

        end

        function setDefaultGUIProperties(obj)
           obj.gui.backgroundColor = [0.5 0.5 0.5];
           obj.gui.elementColor = [0.75 0.75 0.75];
           obj.gui.textColor = [0 0 0];
           
           obj.gui.fontSize = 8;
           obj.gui.fontWeight = 'bold';
           obj.gui.fontName = 'Helvetica';
        end

        function dispDebug(obj,formatSpec,varargin)
            %dispDebug print debug messages (log level >= 4)
            %  input
            %    formatSpec: 	string to print using format specifications similar to fprintf
            %    varargin:   	variables according to formatSpec

            obj.displayToConsole('debug',formatSpec,varargin{:});
        end

        function dispInfo(obj,formatSpec,varargin)
            %dispInfo print information console output (log level >= 3)
            %  input
            %    formatSpec: 	string to print using format specifications similar to fprintf
            %    varargin:   	variables according to formatSpec
            obj.displayToConsole('info',formatSpec,varargin{:});
        end

        function dispError(obj,formatSpec,varargin)
            %dispError print errors (forwarded to "error" that will stop the program) (log level >= 1)
            %  input
            %    formatSpec: 	string to print using format specifications
            %                   similar to 'error'
            %    varargin:   	variables according to formatSpec

            try
                obj.displayToConsole('error',formatSpec,varargin{:});
            catch ME
                if obj.isOctave
                    ME.stack = ME.stack(3:end); % Removes the dispError and dispToConsole from the stack
                    error(ME);
                else
                    throwAsCaller(ME);
                end
            end
        end

        function dispWarning(obj,formatSpec,varargin)
            %dispError print warning (forwarded to 'warning') (log level >= 2)
            %  input
            %    formatSpec: 	string to print using format specifications
            %                   similar to 'warning'
            %    varargin:   	variables according to formatSpec
            obj.displayToConsole('warning',formatSpec,varargin{:});
        end

        function dispDeprecationWarning(obj,formatSpec,varargin)
            %dispDeprecationWarning wrapper for deprecation warnings forwarded to displayToConsole
            obj.displayToConsole('dep',formatSpec,varargin{:});
        end

        function obj = writeLogToFile(obj,filename)
            %writeLogToFile writes the log kept in MatRad_Config to file.
            %  Note that the switch keepLog must be enabled for MatRad_Config to store all logging output.

            singleString = '%s: %s\n';
            fID = fopen(filename,'w');
            fprintf(fID,repmat(singleString,1,size(obj.messageLog,1)),obj.messageLog{:});
            fclose(fID);
        end

        function set.logLevel(obj,newLogLevel)
            %%Property set methods for logLevel
            minLevel = 1;
            maxLevel = 5;
            if newLogLevel >= minLevel && newLogLevel <= maxLevel
                obj.logLevel = newLogLevel;
            else
                obj.dispError('Invalid log level. Value must be between %d and %d',minLevel,maxLevel);
            end
        end

        function set.writeLog(obj,writeLog)
            if writeLog
                logFile = [obj.matRadRoot filesep 'matRad.log'];
                obj.logFileHandle = fopen(logFile,'a');
                obj.writeLog = true;
            else
                fclose(obj.logFileHandle);
                obj.writeLog = false;
            end
        end

        function getEnvironment(obj)
            % getEnvironment function to get the software environment
            %   matRad is running on

            obj.isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
            obj.isMatlab = ~obj.isOctave;

            if obj.isOctave
                obj.env = 'OCTAVE';
                obj.envVersion = OCTAVE_VERSION;
            else
                obj.env = 'MATLAB';
                vData = ver(obj.env);
                obj.envVersion = vData.Version;

            end
        end

        function pln = getDefaultProperties(obj,pln,fields)
            % Function to load all non-set parameters into pln struct
            standardFields = {'propDoseCalc','propOpt','propStf'};

            % Check if only one argument was given
            if ~iscell(fields)
                fields = cellstr(fields);
            end

            for i = 1:length(fields)
                currField = fields{i};

                if ismember(currField,standardFields)
                    % Get defaults for standard fields that can easily be read from set default values
                    if ~isfield(pln,currField)
                        pln.(currField) = struct();
                    end

                    fnames = fieldnames(obj.(currField));
                    for f = 1:length(fnames)
                        if ~isempty(strfind(fnames{f},'default'))
                            cutName = [lower(fnames{f}(8)) fnames{f}(9:end)];
                            if ~isfield(pln.(currField),cutName)
                                pln.(currField).(cutName) = obj.(currField).(fnames{f});
                            end
                        else
                            if ~isfield(pln.(currField),fnames{f})
                                pln.(currField).(fnames{f}) = struct();
                            end
                            subfields = fieldnames(obj.(currField).(fnames{f}));
                            for s = 1:length(subfields)
                                if ~isempty(strfind(subfields{s},'default'))
                                    if length(subfields{s})==8
                                        cutName = [subfields{s}(8)];
                                    else
                                        cutName = [lower(subfields{s}(8)) subfields{s}(9:end)];
                                    end
                                    if ~isfield(pln.(currField).(fnames{f}),cutName)
                                        pln.(currField).(fnames{f}).(cutName) = obj.(currField).(fnames{f}).(subfields{s});
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        function pln = getDefaultClass(obj,pln,propName,configName)
            % load config from pln or from class
            if (isfield(pln,propName) && isstruct(pln.(propName)) && nargin < 4) || (~isfield(pln,propName) && nargin < 4)%if there is no config found
                switch propName
                    case 'propMC'
                        if isfield(pln,'propMC') && strcmp(propName,'propMC') && isfield(pln.propMC,'engine')
                            switch pln.propMC.engine
                                case 'ompMC'
                                    configName = 'matRad_OmpConfig';
                                case 'TOPAS'
                                    configName = 'matRad_TopasConfig';
                                case 'MCsquare'
                                    configName = 'matRad_MCsquareConfig';
                            end
                            pln.propMC = rmfield(pln.propMC,'engine');
                        else
                            if isfield(pln,'radiationMode') && ~isempty(pln.radiationMode)
                                switch pln.radiationMode
                                    case 'photons'
                                        configName = obj.propMC.default_photon_engine;
                                    case 'protons'
                                        configName = obj.propMC.default_proton_engine;
                                    otherwise
                                        configName = obj.propMC.default_carbon_engine;
                                end
                            end
                        end
                    otherwise
                        obj.dispError('Config for ''%s'' not implemented',configName);
                end
            elseif nargin == 4

            elseif nargin < 4 && ~isstruct(pln.(propName))
                % get config name from input field
                configName = class(pln.(propName));
            else
                obj.dispError('Error in default clasee');
            end

            if ~isfield(pln,propName)
                pln.(propName) = struct();
            end
            %Overwrite parameters
            %mc = metaclass(topasConfig); %get metaclass information to check if we can overwrite properties
            if isstruct(pln.(propName))
                % Load configs
                switch configName
                    case 'matRad_OmpConfig'
                        config = matRad_OmpConfig();
                        pln.propMC.engine = 'ompMC';
                    case 'matRad_TopasConfig'
                        config = matRad_TopasConfig();
                        pln.propMC.engine = 'TOPAS';
                    case 'matRad_MCsquareConfig'
                        config = matRad_MCsquareConfig();
                        pln.propMC.engine = 'MCsquare';
                    case 'matRad_HeterogeneityConfig'
                        config = matRad_HeterogeneityConfig;
                end

                props = fieldnames(pln.(propName));
                for fIx = 1:numel(props)
                    fName = props{fIx};
                    if isprop(config,fName)
                        if isstruct(pln.(propName).(fName))
                            SubProps = fieldnames(pln.(propName).(fName));
                            for SubfIx = 1:numel(SubProps)
                                subfName = SubProps{SubfIx};
                                if isfield(config.(fName),subfName)
                                    %We use a try catch block to catch errors when trying
                                    %to overwrite protected/private properties instead of a
                                    %metaclass approach
                                    try
                                        config.(fName).(subfName) = pln.(propName).(fName).(subfName);
                                    catch
                                        obj.dispWarning(['Property ''%s'' for ' configName ' will be omitted due to protected/private access or invalid value.'],fName);
                                    end
                                else
                                    obj.dispWarning(['Unkown property ''%s'' for ' configName ' will be omitted.'],fName);
                                end
                            end
                        else
                            %We use a try catch block to catch errors when trying
                            %to overwrite protected/private properties instead of a
                            %metaclass approach
                            try
                                config.(fName) = pln.(propName).(fName);
                            catch
                                obj.dispWarning(['Property ''%s'' for ' class(config) ' will be omitted due to protected/private access or invalid value.'],fName);
                            end
                        end
                    else
                        obj.dispWarning(['Unkown property ''%s'' for ' class(config) ' will be omitted.'],fName);
                    end
                end

                % Write config to pln
                pln.(propName) = config;
            end

            % Send info to console
            obj.dispInfo(['Class ' class(pln.(propName)) ' has been loaded to pln.' propName '!\n']);

        end

        function configureEnvironment(obj)
            if obj.isOctave
                struct_levels_to_print(0);                  %Disables full printing of struct array fields
                warning("off","Octave:data-file-in-path");  %Disables warning of loading patients from the data folder
            end
        end
    end

    methods(Static)

        function obj = instance()
            %instance creates a singleton instance of MatRad_Config
            %  In MatRad_Config, the constructor is private to make sure only on global instance exists.
            %  Call this static functino to get or create an instance of the matRad configuration class
            persistent uniqueInstance;

            if isempty(uniqueInstance)
                obj = MatRad_Config();
                uniqueInstance = obj;
            else
                obj = uniqueInstance;
            end
        end

        function obj = loadobj(sobj)
            % Overload the loadobj function to allow downward compatibility
            % with workspaces which where saved as an older version of this class

            function basic_struct = mergeStructs(basic_struct, changed_struct)
                % nested function for merging the properties of the loaded
                % obj into a new obj.
                % Merges two structs, including nestes structs, by overwriting
                % the properties of basic_struct with the changed properties in changed_struct
                fields = fieldnames(basic_struct);
                for k = 1:length(fields)
                    disp(fields{k});
                    if(isfield(changed_struct, fields{k}))
                        if isstruct(changed_struct.(fields{k})) && isstruct(basic_struct.(fields{i}))
                            basic_struct.(fields{k}) = mergeStructs(basic_struct.(fields{k}), changed_struct.(fields{i}));
                        else
                            basic_struct.(fields{k}) = changed_struct.(fields{k});
                        end
                    end
                end
            end

            % If the saved object is loaded as a struct there was a problem
            % with the generic loading process most likly a version-conflict
            % regarding the structs, in order to fix this, do a custom
            % loading process including recursivly copying the conflicting structs
            if isstruct(sobj)
                warning('The  loaded object differs from the current MatRad_Config class, resuming the loading process with the overloaded loadobj function!');
                obj = MatRad_Config();
                % Use a metaclass object to get the properties because
                % Octave <= 5.2 doesn't have a properties function
                props = {metaclass(obj).PropertyList.Name};
                % Throw warning if the version differs and remove the
                % matRad_version field from the loaded struct, in order to
                % not overwrite the version later
                if (isfield(sobj, 'matRad_version') && ~(strcmp(obj.matRad_version, sobj.matRad_version)))
                    warning('MatRad version or git Branch of the loaded object differs from the curret version!');
                    sobj = rmfield(sobj, 'matRad_version');
                end
                % Itterate over the properties of the newly created MatRad_Config object
                for i = 1:length(props)
                    % check if the field exists in the loaded object
                    if(isfield(sobj,props{i}))
                        objField = obj.(props{i});
                        sobjField = sobj.(props{i});
                        % If field from loaded object and from the newly
                        % created object are equal skip it, else copy the
                        % value of the loaded object and if it's a struct
                        % check it's field recursively
                        if ~(isequal(sobjField, objField))
                            if (isstruct(sobjField) && isstruct(objField))
                                retStruct = mergeStructs(objField,sobjField);
                                obj.(props{i}) = retStruct;
                            else
                                obj.(props{i}) = sobjField;
                            end
                        end
                    end
                end
            else
                obj = sobj;
            end
        end


    end
end

