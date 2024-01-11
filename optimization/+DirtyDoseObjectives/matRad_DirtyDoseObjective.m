classdef (Abstract) matRad_DirtyDoseObjective < matRad_DoseOptimizationFunction
% matRad_DirtyDoseObjective: Interface for optimization objectives
%   This abstract base class provides the structure of optimization
%   objectives like mean dirty dose, squared deviation dirty dose, EUD dirty dose, dose-volume etc.
%   Implementations can be found in the DirtyDoseObjectives package
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
        penalty                 %Optimization penalty
        
    end
       
    methods (Static)
        function rob = availableRobustness()
            rob = {'none','STOCH','PROB','VWWC','VWWC_INV','COWC','OWC'}; %By default, no robustness is available
        end 
    end
    
    %These should be abstract methods, however Octave can't parse them. As soon
    %as Octave is able to do this, they should be made abstract again
    methods %(Abstract)
       
        %returns the objective function value for the given dose vector. Needs to be implemented in sub-classes.
        function fDirtyDose       = computeDirtyDoseObjectiveFunction(obj,dirtyDose)
            error('Function needs to be implemented!');
        end
        
        
        %returns the dose-gradient for the given dose vector. Needs to be implemented in sub-classes.
        function fDirtyDoseGrad   = computeDirtyDoseObjectiveGradient(obj,dirtyDose)
            error('Function needs to be implemented!');
        end

    end
    
    methods (Access = public)
       
        % constructor of matRad_DoseObjective
        function obj = matRad_DirtyDoseObjective(varargin)
            %default initialization from struct (parameters & penalty)
            obj@matRad_DoseOptimizationFunction(varargin{:});
        end
        
        %Overloads the struct function to add Objective related information
        %to output struct
        function s = struct(obj)
            s = struct@matRad_DoseOptimizationFunction(obj);
            s.penalty = obj.penalty;
            
        end
        
    end   
end

