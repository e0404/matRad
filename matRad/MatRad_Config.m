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
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
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

        defaults;

        %Disable GUI
        disableGUI = false;

        devMode = false;
        eduMode = false;
        
        gui;

        %User folders
        userfolders; %Cell array of user folders containing machines, patients, hluts. Default contains the userdata folder in the matRad root directory
    end
    
    %Deprecated properties referencing a newer one
    properties (Dependent,SetAccess = private)
        %Default Properties
        propDoseCalc;
        propOpt;
        propStf;
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

    properties (SetAccess = private, Dependent)
        matRadSrcRoot;      %Path to matRadSrcRoot ("matRad" subfolder of matRadRoot)
        primaryUserFolder;  %Points to the first entry in userfolders
        exampleFolder;      %Contains examples  
        thirdPartyFolder;   %Contains third party tools/libraries used in matRad
    end

    methods (Access = private)
        function obj = MatRad_Config()
            %MatRad_Config Constructs an instance of this class.
            %  The configuration is implemented as a singleton and used globally
            %  Therefore its constructor is private
            %  For instantiation, use the static MatRad_Config.instance();
            
            %Set Path
            if isdeployed
                obj.matRadRoot = [ctfroot filesep 'matRad'];

                if ispc
                    userdir= getenv('USERPROFILE');
                else 
                    userdir= getenv('HOME');
                end

                userfolderInHomeDir = [userdir filesep 'matRad'];               

                obj.userfolders = {userfolderInHomeDir};
            else
                obj.matRadRoot = fileparts(fileparts(mfilename('fullpath')));
                addpath(genpath(obj.matRadSrcRoot));
                addpath(obj.exampleFolder);
                addpath(genpath(obj.thirdPartyFolder));
                obj.userfolders = {[obj.matRadRoot filesep 'userdata' filesep]};
            end           
            
            %set version
            obj.getEnvironment();
            obj.matRad_version = matRad_version(obj.matRadRoot);

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
                        warning('matRad:Deprecated',forwardArgs{:});
                    end
                case{'warning'}
                    if obj.logLevel >= 2
                        warning('matRad:Warning',forwardArgs{:});
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
            
            %Default Steering/Geometry Properties
            obj.defaults.propStf.generator = {'PhotonIMRT','ParticleIMPT','SimpleBrachy'};
            obj.defaults.propStf.longitudinalSpotSpacing = 2;
            obj.defaults.propStf.addMargin = true; %expand target for beamlet finding
            obj.defaults.propStf.bixelWidth = 5;
          
            %Dose Calculation Options
            obj.defaults.propDoseCalc.engine = {'SVDPB','HongPB'}; %Names for default engines used when no other is given
            obj.defaults.propDoseCalc.doseGrid.resolution = struct('x',3,'y',3,'z',3); %[mm]
            obj.defaults.propDoseCalc.dosimetricLateralCutOff = 0.995; %[rel.]
            obj.defaults.propDoseCalc.geometricLateralCutOff = 50; %[mm]
            obj.defaults.propDoseCalc.kernelCutOff = Inf; %[mm]
            obj.defaults.propDoseCalc.ssdDensityThreshold = 0.05; %[rel.]
            obj.defaults.propDoseCalc.useGivenEqDensityCube = false; %Use the given density cube ct.cube and omit conversion from cubeHU.
            obj.defaults.propDoseCalc.ignoreOutsideDensities = true; %Ignore densities outside of cst contours
            obj.defaults.propDoseCalc.useCustomPrimaryPhotonFluence = false; %Use a custom primary photon fluence
            obj.defaults.propDoseCalc.calcLET = true; %calculate LETs for particles
            obj.defaults.propDoseCalc.selectVoxelsInScenarios = 'all';
            obj.defaults.propDoseCalc.airOffsetCorrection = true;
            % default properties for fine sampling calculation
            obj.defaults.propDoseCalc.fineSampling.sigmaSub = 1;
            obj.defaults.propDoseCalc.fineSampling.N = 2;
            obj.defaults.propDoseCalc.fineSampling.method = 'fitCircle';
            %Monte Carlo options
            obj.defaults.propDoseCalc.numHistoriesPerBeamlet = 2e4;
            obj.defaults.propDoseCalc.numHistoriesDirect = 1e6;
            obj.defaults.propDoseCalc.outputMCvariance = true;
                      
            %Optimization Options
            obj.defaults.propOpt.optimizer = 'IPOPT';
            obj.defaults.propOpt.maxIter = 500;
            obj.defaults.propOpt.runDAO = 0;
            obj.defaults.propOpt.clearUnusedVoxels = false;

            %Sequencing Options
            obj.defaults.propSeq.sequencer = 'siochi';
            


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

            obj.logLevel   = 1; %Omit output except errors
            
            %Default Steering/Geometry Properties
            obj.defaults.propStf.longitudinalSpotSpacing = 20;
            obj.defaults.propStf.bixelWidth = 20;
            
            %Dose Calculation Options
            obj.defaults.propDoseCalc.doseGrid.resolution = struct('x',5,'y',6,'z',7); %[mm]
            obj.defaults.propDoseCalc.geometricLateralCutOff = 20;
            obj.defaults.propDoseCalc.dosimetricLateralCutOff = 0.8;
            obj.defaults.propDoseCalc.kernelCutOff = 20; %[mm]

            %Monte Carlo options
            obj.defaults.propDoseCalc.numHistoriesPerBeamlet = 100;
            obj.defaults.propDoseCalc.numHistoriesDirect = 100;
            
            %Optimization Options
            obj.defaults.propOpt.maxIter = 10;

            obj.defaults.samplingScenarios = 2;

            obj.disableGUI = true;

            obj.devMode = true;
            obj.eduMode = false;
        end  
        
        %%for edu mode
        function setDefaultPropertiesForEduMode(obj)
            obj.setDefaultProperties();

            obj.logLevel = 2;
            
            %Default Steering/Geometry Properties
            obj.defaults.propStf.longitudinalSpotSpacing = 3;
            
            %Dose calculation options
            obj.defaults.propDoseCalc.resolution = struct('x',4,'y',4,'z',4); %[mm]
            obj.defaults.propDoseCalc.lateralCutOff = 0.975; %[rel.]
            
            %Optimization Options
            obj.defaults.propOpt.maxIter = 500;
                       
            obj.disableGUI = false;
            
            obj.devMode = false;
            obj.eduMode = true;
        end

        function setDefaultGUIProperties(obj)
            %Detect current theme
            light = false;
            try
                if ispc
                    light = logical(winqueryreg('HKEY_CURRENT_USER','Software\\Microsoft\\Windows\\CurrentVersion\\Themes\\Personalize','AppsUseLightTheme'));
                elseif ismac
                    out = system('defaults read -g AppleInterfaceStyle');
                    if ~strcmp(out,'Dark')
                        light = true;
                    end
                else
                    out = system('gsettings get org.gnome.desktop.interface color-scheme');
                    if strcmp(out,'prefer-light')
                        light = true;
                    end
                end
            catch
                light = false;
            end
            
            if light
                theme = matRad_ThemeLight();
            else
                theme = matRad_ThemeDark();
            end

            obj.gui = struct(theme);
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

        function set.userfolders(obj,userfolders)
            oldFolders = obj.userfolders;
                     
            %Check if folders need to be created
            for f = 1:numel(userfolders)
                if ~isfolder(userfolders{f})
                    [status, msg] = mkdir(userfolders{f});
                    if status == 0
                        obj.dispWarning('Userfolder %s not added beacuse it could not be created: %s',userfolders{f},msg);
                    else
                        subfolders = {'hluts','machines','patients','scripts'};                    
                        [status,msgs] = cellfun(@(sub) mkdir([userfolders{f} filesep sub]),subfolders,'UniformOutput',false);
                        if any(cell2mat(status) ~= 1)
                            obj.dispWarning('Problem when creating subfolder in Userfolder %s!',userfolders{f})
                        end
                    end
                end
            end

            %We do this to verify folders
            nonWorkingFolders = cellfun(@isempty,userfolders);
            userfolders(nonWorkingFolders) = [];

            allNewFolders = cellfun(@dir, userfolders,'UniformOutput',false);
            if isempty(allNewFolders)
                obj.dispWarning('No user folders specified. Defaulting to userdata folder in matRad root directory.');
                if ~isdeployed
                    allNewFolders = {[fileparts(mfilename('fullpath')) filesep 'userdata' filesep]}; %We don't access obj.matRadRoot here because of Matlab's weird behavior with properties
                else
                    allNewFolders = {[ctfroot filesep 'userdata' filesep]}; %We don't access obj.matRadRoot here because of Matlab's weird behavior with properties
                end
            end           

            cleanedNewFolders = cellfun(@(x) x(1).folder,allNewFolders,'UniformOutput',false);
            
            % Identify newly added folder paths and add them to path
            if ~isdeployed
                if ~isempty(oldFolders) %if statement for octave compatibility
                    addedFolders = setdiff(cleanedNewFolders, oldFolders);
                else
                    addedFolders = cleanedNewFolders;
                end
                addedFolders = cellfun(@genpath,addedFolders,'UniformOutput',false);
                addedFolders = strjoin(addedFolders,pathsep);
                addpath(addedFolders);
            end

            % Identify removed folder paths
            if ~isempty(oldFolders) %if statement for octave compatibility
                removedFolders = setdiff(oldFolders, cleanedNewFolders);
                removedFolders = cellfun(@genpath,removedFolders,'UniformOutput',false);
                removedFolders = strjoin(removedFolders,pathsep);
                if ~isdeployed
                    rmpath(removedFolders);
                end
            end
            
            obj.userfolders = cleanedNewFolders;
        end

        function srcRoot = get.matRadSrcRoot(obj)
            srcRoot = [obj.matRadRoot filesep 'matRad' filesep];
        end

        function primaryUserFolder = get.primaryUserFolder(obj)
            primaryUserFolder = obj.userfolders{1};
        end

        function exampleFolder = get.exampleFolder(obj)
            exampleFolder = [obj.matRadRoot filesep 'examples' filesep];            
        end

        function thirdPartyFolder = get.thirdPartyFolder(obj)
            thirdPartyFolder = [obj.matRadRoot filesep 'thirdParty' filesep];
        end

        function propDoseCalc = get.propDoseCalc(obj)
            obj.dispWarning('Property ''propDoseCalc'' is deprecated. Use ''defaults.propDoseCalc'' instead!');
            
            fNames = fieldnames(obj.defaults.propDoseCalc);            
            for i = 1:numel(fNames)
                fNewName = ['default' upper(fNames{i}(1)) fNames{i}(2:end)];
                propDoseCalc.(fNewName) = obj.defaults.propDoseCalc.(fNames{i});
            end
        end

        function propStf = get.propStf(obj)
            obj.dispWarning('Property ''propStf'' is deprecated. Use ''defaults.propStf'' instead!');
            
            fNames = fieldnames(obj.defaults.propStf);            
            for i = 1:numel(fNames)
                fNewName = ['default' upper(fNames{i}(1)) fNames{i}(2:end)];
                propStf.(fNewName) = obj.defaults.propStf.(fNames{i});
            end
        end

        function propOpt = get.propOpt(obj)
            obj.dispWarning('Property ''propOpt'' is deprecated. Use ''defaults.propStf'' instead!');
            
            fNames = fieldnames(obj.defaults.propOpt);            
            for i = 1:numel(fNames)
                fNewName = ['default' upper(fNames{i}(1)) fNames{i}(2:end)];
                propOpt.(fNewName) = obj.defaults.propOpt.(fNames{i});
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
                        pln.(currField) = obj.defaults.(currField);
                    else
                        pln.(currField) = matRad_recursiveFieldAssignment(pln.(currField),obj.defaults.(currField),false);
                    end
                end
            end
        end

        function configureEnvironment(obj)
            if obj.isOctave
                struct_levels_to_print(0);                  %Disables full printing of struct array fields
                warning("off","Octave:data-file-in-path");  %Disables warning of loading patients from the data folder
            end
        end
    end

    %methods (Access = private)
    %    function renameFields(obj)
    %end

    methods (Static)

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
                meta = metaclass(obj);
                props = {meta.PropertyList.Name};
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

