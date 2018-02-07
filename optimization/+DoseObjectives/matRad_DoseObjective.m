classdef (Abstract) matRad_DoseObjective
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_DoseObjective: Interface for optimization objectives 
%   This abstract base class provides the structure of optimization
%   objectives like mean dose, squared deviation, EUD, dose-volume etc.
%   Implementations can be found in the DoseObjectives package
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
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
        parameters          %Cell array of parameter values       
        penalty             %Optimization penalty
    end
    
    methods (Abstract)
        %returns the objective function value for the given dose vector. Needs to be implemented in sub-classes.
        fDose       = computeDoseObjectiveFunction(obj,dose) 
        
        %returns the dose-gradient for the given dose vector. Needs to be implemented in sub-classes.
        fDoseGrad   = computeDoseObjectiveGradient(obj,dose) 
    end        
    
    %Helper methods
    methods (Access = public)
        function doseParams = getDoseParameters(obj)
            %Get only the dose related parameters.
            ix = cellfun(@(c) isequal('dose',c),obj.parameterTypes);
            doseParams = [obj.parameters{ix}];
        end
        
        function obj = setDoseParameters(obj,doseParams)
            %Set only the dose related parameters.
            ix = cellfun(@(c) isequal('dose',c),obj.parameterTypes);
            obj.parameters(ix) = num2cell(doseParams);

        end
    end
end

