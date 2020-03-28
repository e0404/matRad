classdef MatRad_Config < handle
    %MatRad_Config MatRad Configuration
    
    properties
        logLevel = 3;
        propDoseCalc;        
        propOpt;
        propMC;
        propStf;
        keepLog = false;        
        %logToFile = [];
    end
    
    properties (SetAccess = private)
        messageLog = {};
    end
    
    properties (Constant)
        matRadRoot = fileparts(mfilename('fullpath'));
    end
    
    methods (Access = private)
        function obj = MatRad_Config()
            %MatRad_Config Construct an instance of this class

            obj.propStf.defaultLongitudinalSpotSpacing = 3;
            
            obj.propDoseCalc.defaultResolution = struct('x',3,'y',3,'z',3); %[mm]
            obj.propDoseCalc.defaultLateralCutOff = 0.995; %[rel.]
            obj.propDoseCalc.defaultGeometricCutOff = 50; %[mm]
            obj.propDoseCalc.ssdDensityThreshold = 0.05; %[rel.]
            
            obj.propOpt.defaultMaxIter = 500;
            
            obj.propMC.ompMC_defaultHistories = 1e6;
            obj.propMC.MCsquare_defaultHistories = 1e6;
            obj.propMC.direct_defaultHistories = 2e4;
            

        end
                
    end
    
    methods
        function dispDebug(obj,formatSpec,varargin)
            obj.displayToConsole('debug',formatSpec,varargin{:});
        end
        
        function dispInfo(obj,formatSpec,varargin)
            obj.displayToConsole('info',formatSpec,varargin{:});
        end
        
        function dispError(obj,formatSpec,varargin)
            obj.displayToConsole('error',formatSpec,varargin{:});
        end
        
        function dispWarning(obj,formatSpec,varargin)
            obj.displayToConsole('warning',formatSpec,varargin{:});
        end
        
        function displayToConsole(obj,type,formatSpec,varargin)
            if nargin < 4
                forwardArgs = {formatSpec};
            else
                forwardArgs = {formatSpec,varargin{:}};
            end
            
            if obj.keepLog
                obj.messageLog{end+1,1} = upper(type);
                obj.messageLog{end,2} = sprintf(forwardArgs{:});
            end
            
            switch type
                case {'error'}
                    if obj.logLevel >= 1                            
                        error(forwardArgs{:});                        
                    end
                case{'warning'}
                    if obj.logLevel >= 2
                        warning(forwardArgs{:});
                    end
                case{'info'}                    
                    if obj.logLevel >= 3
                        fprintf(forwardArgs{:});                    
                    end
                case{'debug'}                    
                    if obj.logLevel >= 4
                        forwardArgs{1} = ['DEBUG: ' forwardArgs{1}];
                        fprintf(forwardArgs{:});
                    end
                otherwise
                    error('Log type %s not defined!',type);
            end
        end
        
        function obj = writeLogToFile(obj,filename)
            singleString = '%s: %s\n';
            fID = fopen(filename,'w');
            fprintf(fID,repmat(singleString,1,size(obj.messageLog,1)),obj.messageLog{:});
            fclose(fID);
        end
        
        %%Property set methods for checks
        function set.logLevel(obj,newLogLevel)
            minLevel = 1;
            maxLevel = 4;
            if newLogLevel >= minLevel && newLogLevel <= maxLevel
                obj.logLevel = newLogLevel;
            else
                obj.dispError('Invalid log level. Value must be between %d and %d',minLevel,maxLevel);
            end
        end
    end
    
    methods(Static)
        function obj = instance()
            persistent uniqueInstance;
            
            if isempty(uniqueInstance)
                obj = MatRad_Config();
                uniqueInstance = obj;
            else
                obj = uniqueInstance;
            end
        end
    end
end

