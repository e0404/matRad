classdef matRad_OptimizationProblemDAO < matRad_OptimizationProblem
% handle class to keep state easily
%
% References
%   -
%
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
    
    % properties (Access = private)
    %    currentDose
    %    currentWeights
    %end
    properties
        apertureInfo
    end    
    
    methods (Static)
        %In External Files
        updatedInfo = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfoVect);
        
        [apertureInfoVec, mappingMx, limMx] = matRad_daoApertureInfo2Vec(apertureInfo);
    end
    
    methods
        function obj = matRad_OptimizationProblemDAO(backProjection,apertureInfo)
            obj = obj@matRad_OptimizationProblem(backProjection);
            obj.apertureInfo = apertureInfo;
        end       
        
        function lb = lowerBounds(obj,w)
            lb = obj.apertureInfo.limMx(:,1);  % Lower bound on the variables. 
        end
        
        function ub = upperBounds(obj,w)            
            ub = obj.apertureInfo.limMx(:,2);
        end
    end
end

