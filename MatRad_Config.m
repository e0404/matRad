classdef MatRad_Config < handle
    %MatRad_Config MatRad Configuration
    
    properties
        logLevel = 1;
        propDoseCalc;
        messageLog = {};
        keepLog = false;
        %logToFile = [];
    end
    
    methods (Access = private)
        function obj = MatRad_Config()
            %MatRad_Config Construct an instance of this class
            
            obj.propDoseCalc.defaultResolution = struct('x',3,'y',3,'z',3); %[mm]
            obj.propDoseCalc.defaultLateralCutOff = 0.995; %[rel.]
            obj.propDoseCalc.defaultGeometricCutOff = 50; %[mm]
            obj.propDoseCalc.ssdDensityThreshold = 0.05; %[rel.]
        end
                
    end
    
    methods
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
                    forwardArgs{1} = [forwardArgs{1} '\n'];
                    if obj.logLevel >= 3
                        fprintf(forwardArgs{:});                    
                    end
                case{'debug'}
                    %forwardArgs{1}(end+1) = '\n';                    
                    if obj.logLevel >= 4
                        forwardArgs{1} = ['DEBUG: ' forwardArgs{1} '\n'];
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

