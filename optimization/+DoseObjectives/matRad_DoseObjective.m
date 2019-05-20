classdef (Abstract) matRad_DoseObjective < matRad_DoseOptimizationFunction
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
    
    properties (Abstract, Access = public)
        penalty             %Optimization penalty
    end
    
    %These should be abstract methods, however Octave can't parse them. As soon
    %as Octave is able to do this, they should be made abstract again
    methods %(Abstract)
       
        %returns the objective function value for the given dose vector. Needs to be implemented in sub-classes.
        function fDose       = computeDoseObjectiveFunction(obj,dose)
            error('Function needs to be implemented!');
        end
        
        
        %returns the dose-gradient for the given dose vector. Needs to be implemented in sub-classes.
        function fDoseGrad   = computeDoseObjectiveGradient(obj,dose)
            error('Function needs to be implemented!');
        end
         
    end
    
    methods (Access = public)
       
       % default constructor of matRad_DoseObjective
        function obj = matRad_DoseObjective(inputArg)
            % superclass constructor is already called when this is line is reached
            % additional matRad_DoseObjective constructor specific code goes here
        end
        
        %Overloads the struct function to add Objective related information
        %to output struct
        function s = struct(obj)
            s = struct@matRad_DoseOptimizationFunction(obj);
            s.penalty = obj.penalty;
        end
    end 
end

