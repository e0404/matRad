classdef (Abstract) matRad_PhantomVOIVolume < handle
% matRad_PhantomVOIVolume: Interface for VOI Volumes
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

    methods
        function obj = matRad_PhantomVOIVolume(name,type,p)
        %p is the input parser used in the child classes to check for additional variables
        
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
            nxIdx = size(cst,1)+1;
            cst{nxIdx,1}                = nxIdx-1;
            cst{nxIdx,2}                = obj.name;
            cst{nxIdx,3}                = obj.type;
            cst{nxIdx,5}.TissueClass    = obj.TissueClass;
            cst{nxIdx,5}.alphaX         = obj.alphaX;
            cst{nxIdx,5}.betaX          = obj.betaX;
            cst{nxIdx,5}.Priority       = nxIdx;
            cst{nxIdx,5}.Visible        = obj.Visible;

            if nxIdx <= size(obj.colors,1)
                obj.visibleColor = obj.colors(nxIdx,:);
            end
            cst{nxIdx,5}.visibleColor   = obj.visibleColor;

            if ~iscell(obj.objectives) %should be redundant
                DoseObjectives = {obj.objectives};  
            else
                DoseObjectives = obj.objectives;
            end
            for i = 1:numel(DoseObjectives)
                cst{nxIdx,6} {i}= DoseObjectives{i}; 
            end
        end
    end
end