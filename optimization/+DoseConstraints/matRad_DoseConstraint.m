classdef (Abstract) matRad_DoseConstraint
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_DoseConstraint: Interface for optimization constraints.
%   This abstract base class provides the interface of constraints for
%   non-linear optimization.
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
    end
        
    methods (Abstract)
        %returns the constraint function(s) value(s) for a given dose
        %vector. Needs to be implemented in sub-classes.
        cDose        = computeDoseConstraintFunction(obj,dose)
        
        %return the (dose-dependent) constraint function jacobian for a
        %given dose vector. Needs to be implemented in sub-classes.
        cDoseJacob   = computeDoseConstraintJacobian(obj,dose)
        
        %Returns upper bound(s) / max value(s) for constraint function(s)
        %Needs to be implemented in sub-classes.
        cu           = upperBounds(obj,n)
        
        %Returns lower bound(s) / min value(s) for constraint function(s)
        %Needs to be implemented in sub-classes.
        cl           = lowerBounds(obj,n)                
    end

    methods (Access = public)
        function jStruct = getDoseConstraintJacobianStructure(obj,n)
        %return the structure of the (dose-dependent) constraint function 
        %jacobian for a given length n of the dose vector. Returns a
        %default of a jStruct
            jStruct = ones(n,1);
        end
  
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

