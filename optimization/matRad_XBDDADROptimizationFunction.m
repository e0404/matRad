classdef (Abstract) matRad_XBDDADROptimizationFunction
% matRad_DoseOptimizationFunction. Superclass for objectives and constraints
% This is the superclass for all objectives and constraints to enable easy 
% one-line identification.
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
    
    properties (Abstract, Constant)
        name                %Display name of the Objective. Needs to be implemented in sub-classes.
        parameterNames      %Cell array of Display names of the parameters. Needs to be implemented in sub-classes.
        parameterTypes      %Cell array of parameter types. Valid types are 'dose', 'numeric', or a cell list of string options. Needs to be implemented in sub-classes.
    end
    
    properties (Abstract, Access = public)
        parameters
    end
    
    methods
        function obj = matRad_XBDDADROptimizationFunction(dataStruct)
            if nargin > 0 && ~isempty(dataStruct) && isstruct(dataStruct)
                obj = assignCommonPropertiesFromStruct(obj,dataStruct);
            end
        end
        
        % Overload the struct function to return a struct with general
        % the objective / constraint
        function s = struct(obj)
            s.className = class(obj);
            s.parameters = obj.parameters;
        end
    end
    
    % Helper methods
    methods (Access = public)
        function doseParams = getXBDDADRParameters(obj)
            % get only the dose related parameters.
            ix = cellfun(@(c) isequal('XBD(DADR)',c),obj.parameterTypes);
            doseParams = [obj.parameters{ix}];
        end
        
        function obj = setXBDDADRParameters(obj,doseParams)
            % set only the dose related parameters.
            ix = cellfun(@(c) isequal('XBD(DADR)',c),obj.parameterTypes);
            obj.parameters(ix) = num2cell(doseParams);
            
        end                
    end
    
    methods (Access = private)
        function obj = assignCommonPropertiesFromStruct(obj,s)
            for fn = fieldnames(s)'    %enumerat fields
                try
                    obj.(fn{1}) = s.(fn{1});   %and copy
                catch
                    continue;
                    %Do Nothing here
                    %warning('Could not copy field %s', fn{1});
                end
            end
        end
    end
        
    
    methods (Static)
        % creates an optimization function from a struct
        function obj = createInstanceFromStruct(s)            
            try
                %Check vor old version of cst objectives / cosntraints and
                %convert if necessary
                if isfield(s,'type')
                    s = matRad_DADROptimizationFunction.convertOldOptimizationStruct(s);
                end
                
                %Create objective / constraint from class name
                obj = eval([s.className '(s)']);       
                
                env = matRad_getEnvironment();
                
                %Workaround for Octave which has a problem assigning
                %properties in superclass
                if strcmp(env,'OCTAVE')
                    obj = assignCommonPropertiesFromStruct(obj,s);
                end
                
            catch ME
                error(['Could not instantiate Optimization Function: ' ME.message]);
            end
        end                
    end
end