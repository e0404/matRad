classdef (Abstract) matRad_VOIVolume < handle
% matRad_VOIVolume: Interface for VOI Volumes
%   This abstract base class provides the structure of VOI Volumes.
%   So far implemented: Box and spherical objectives   
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team.
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
        idx;
        name;
        type;
        TissueClass = 1;
        alphaX = 0.1000;
        betaX = 0.0500;
        Priority = 1;
        Visible = 1;
        visibleColor = [0 0 0];
        HU = 0;
        offset = [0,0,0]; %center of objective
        objectives = {};
        colors = [[1,0,0];[0,1,0];[0,0,1];[1,1,0];[1,0,1];[0,1,1];[1,1,1]];
    end

    methods (Static, Access = private)

        function oldValue = getOrIncrementCount(increment) %used to automatically index the objectives
        % Private function to manage the counter
            persistent VALUE
            if isempty(VALUE)
                VALUE = 0;
            end
            oldValue = VALUE;
            if nargin > 0
                VALUE = VALUE + increment;
            end
        end 
    end

     methods (Static)
        function value = getInstanceCount()
        % Public access to the counter cannot increment it
            value = cldef.getOrIncrementCount();
        end
    end

    methods
        function obj = matRad_VOIVolume(name,type,p)
        %p is the input parser used in the child classes to check for additional variables
        % Increment the counter in the constructor
            matRad_VOIVolume.getOrIncrementCount(1);
            obj.idx = matRad_VOIVolume.getOrIncrementCount();
            if obj.idx <= size(obj.colors,1)
                obj.visibleColor = obj.colors(obj.idx,:);
            else
                obj.visibleColor = [1 1 1];
            end
            obj.Priority = obj.idx;
            obj.name = name;
            obj.type = type;
            obj.offset = p.Results.offset;
            obj.HU = p.Results.HU;


            %idea is that DoseObjectiveFunction can be either a single objective or an 
            %array of objectives. If it is a single objective store it as a cell array 
            if iscell(p.Results.objectives)
                obj.objectives = p.Results.objectives;
            else 
                obj.objectives = {p.Results.objectives};
            end
            %}
        end

        function cst = initializeParameters(obj,cst)
            %initialize entry for this VOI in cst
            cst{obj.idx,1}                = obj.idx-1;
            cst{obj.idx,2}                = obj.name;
            cst{obj.idx,3}                = obj.type;
            cst{obj.idx,5}.TissueClass    = obj.TissueClass;
            cst{obj.idx,5}.alphaX         = obj.alphaX;
            cst{obj.idx,5}.betaX          = obj.betaX;
            cst{obj.idx,5}.Priority       = obj.Priority;
            cst{obj.idx,5}.Visible        = obj.Visible;
            cst{obj.idx,5}.visibleColor   = obj.visibleColor;

            if ~iscell(obj.objectives) %should be redundant
                DoseObjectives = {obj.objectives};  
            else
                DoseObjectives = obj.objectives;
            end
            for i = 1:numel(DoseObjectives)
                cst{obj.idx,6} {i}= DoseObjectives{i}; 
            end
        end
    end
end