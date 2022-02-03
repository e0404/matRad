classdef (Abstract) matRad_LETObjective < matRad_LETOptimizationFunction
% matRad_LETObjective: Interface for optimization objectives
%   This abstract base class provides the structure of optimization
%   objectives like mean LET, squared deviation, EULET, LET-volume etc.
%   Implementations can be found in the LETObjectives package
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team.
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
       
        %returns the objective function value for the given LET vector. Needs to be implemented in sub-classes.
        function fLET       = computeLETObjectiveFunction(obj,LET)
            error('Function needs to be implemented!');
        end
        
        
        %returns the LET-gradient for the given LET vector. Needs to be implemented in sub-classes.
        function fLETGrad   = computeLETObjectiveGradient(obj,LET)
            error('Function needs to be implemented!');
        end
         
    end
    
    methods (Access = public)
       
        % constructor of matRad_LETObjective
        function obj = matRad_LETObjective(varargin)
            %default initialization from struct (parameters & penalty)
            obj@matRad_LETOptimizationFunction(varargin{:});
        end
        
        %Overloads the struct function to add Objective related information
        %to output struct
        function s = struct(obj)
            s = struct@matRad_LETOptimizationFunction(obj);
            s.penalty = obj.penalty;
        end
    end 
end

