classdef MatRad_Config < handle
    %MatRad_Config MatRad Configuration class
    % This class is used globally through Matlab to handle default values and 
    % logging and is declared as global matRad_cfg.
    % Usage:
    %   global matRad_cfg; matRad_cfg = MatRadConfig.instance();    
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
        logLevel = 3;
        propDoseCalc;        
        propOpt;
        propMC;
        propStf;
        keepLog = false; 
        
        disableGUI = false;
    end
    
    properties (SetAccess = private)
        messageLog = {};
    end
    
    properties (Constant)
        matRadRoot = fileparts(mfilename('fullpath'));
    end
    
    methods (Access = private)
        function obj = MatRad_Config()
            %MatRad_Config Constructs an instance of this class.             
            %  The configuration is implemented as a singleton and used globally
            %  Therefore its constructor is private
            %  For instantiation, use the static MatRad_Config.instance();           
            obj.setDefaultProperties();           
        end
                
    end
    
    methods
        function setDefaultProperties(obj)
            obj.propStf.defaultLongitudinalSpotSpacing = 3;
            obj.propStf.defaultAddMargin = true; %expand target for beamlet finding
            
            obj.propDoseCalc.defaultResolution = struct('x',3,'y',3,'z',3); %[mm]
            obj.propDoseCalc.defaultLateralCutOff = 0.995; %[rel.]
            obj.propDoseCalc.defaultGeometricCutOff = 50; %[mm]
            obj.propDoseCalc.ssdDensityThreshold = 0.05; %[rel.]            
            obj.propOpt.defaultMaxIter = 500;           
            obj.propMC.ompMC_defaultHistories = 1e6;
            obj.propMC.ompMC_defaultOutputVariance = false;
            obj.propMC.MCsquare_defaultHistories = 1e6;
            obj.propMC.direct_defaultHistories = 2e4;
            obj.disableGUI = false;
        end
  
        %%FOr testing
        function setDefaultPropertiesForTesting(obj)
            obj.logLevel   = 1;
            obj.propStf.defaultLongitudinalSpotSpacing = 20;
            obj.propStf.defaultAddMargin = true; %expand target for beamlet finding
            obj.propDoseCalc.defaultResolution = struct('x',5,'y',6,'z',7); %[mm]
            obj.propDoseCalc.defaultGeometricCutOff = 20;
            obj.propDoseCalc.defaultLateralCutOff = 0.8;
            obj.propOpt.defaultMaxIter = 10;
            obj.propMC.ompMC_defaultHistories = 100;
            obj.propMC.MCsquare_defaultHistories = 100;
            obj.propMC.direct_defaultHistories = 100;
            obj.disableGUI = true;
        end  
      
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

