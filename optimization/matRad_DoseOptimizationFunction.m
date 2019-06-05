classdef (Abstract) matRad_DoseOptimizationFunction
% matRad_DoseOptimizationFunction This is the superclass for all
% objectives and constraints to enable easy one-line identification
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
        function obj = matRad_DoseOptimizationFunction(dataStruct)
            if nargin == 0
                return;
            end
            if isstruct(dataStruct)
                for fn = fieldnames(dataStruct)'    %enumerat fields
                    try
                        obj.(fn{1}) = dataStruct.(fn{1});   %and copy
                    catch
                        %Do Nothing here
                        %warning('Could not copy field %s', fn{1});
                    end
                end
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
        function doseParams = getDoseParameters(obj)
            % get only the dose related parameters.
            ix = cellfun(@(c) isequal('dose',c),obj.parameterTypes);
            doseParams = [obj.parameters{ix}];
        end
        
        function obj = setDoseParameters(obj,doseParams)
            % set only the dose related parameters.
            ix = cellfun(@(c) isequal('dose',c),obj.parameterTypes);
            obj.parameters(ix) = num2cell(doseParams);
            
        end
    end
    
    methods (Static)
        % creates an optimization function from a struct
        function obj = createInstanceFromStruct(s)
            try
                obj = eval([s.className '(s)']);
            catch
                error(['Input struct / Parameter invalid for creation of optimization function!']);
            end
        end
    end
end

