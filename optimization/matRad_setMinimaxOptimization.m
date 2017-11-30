function [cst,dij,pln] = matRad_setMinimaxOptimization(cst,dij,pln)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad intitialization for minimax or maximin optimization. Adds an
% auxiliary variable and constraint for each minimax/maximin objective. The
% auxiliary variable is the max/min dose that will be optimized, while the 
% auxiliary constraint makes sure that no structure voxels exceed this 
% max/min dose
% 
% call
%   [cst,dij,pln] = matRad_setMinimaxOptimization(cst,dij,pln)
%
% input
%   cst:        matRad cst struct
%   dij:        matRad dij struct
%   pln:        matRad pln struct
%
% output
%   cst:        matRad cst struct
%   dij:        matRad dij struct
%   pln:        matRad pln struct
%
% References
%   Boyd and L. Vandenberghe, Convex Optimization (Cambridge University 
%   Press, Cambridge, UK, 2004
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize auxiliary variable number
auxVarNumber = 0;

% loop over objectives/constraints
for  i = 1:size(cst,1)
    for j = 1:numel(cst{i,6})
        if isequal(cst{i,6}(j).type, 'max dose objective (exact)') || isequal(cst{i,6}(j).type, 'min dose objective (exact)')

            % set auxiliary variable number
            auxVarNumber = auxVarNumber + 1;
            cst{i,6}(j).auxVarNum = dij.totalNumOfBixels + auxVarNumber;

            % add auxiliary constraint
            cst{i,6}(numel(cst{i,6}) + 1) = cst{i,6}(j);
            if isequal(cst{i,6}(j).type, 'max dose objective (exact)')
                cst{i,6}(numel(cst{i,6})).type = 'minimax constraint (exact)';
                cst{i,6}(numel(cst{i,6})).dose = 0;
                cst{i,6}(numel(cst{i,6})).penalty = NaN;
                cst{i,6}(numel(cst{i,6})).auxVarNum = dij.totalNumOfBixels + auxVarNumber;
            elseif isequal(cst{i,6}(j).type, 'min dose objective (exact)')
                cst{i,6}(numel(cst{i,6})).type = 'maximin constraint (exact)';
                cst{i,6}(numel(cst{i,6})).dose = 0;
                cst{i,6}(numel(cst{i,6})).penalty = NaN;
                cst{i,6}(numel(cst{i,6})).auxVarNum = dij.totalNumOfBixels + auxVarNumber;
            end                
        end
    end
end

% add auxiliary variables
dij.physicalDose{1}(:,end+1:end+auxVarNumber) = sparse(size(dij.physicalDose{1},1), auxVarNumber);
dij.totalNumOfBixels = dij.totalNumOfBixels + auxVarNumber;
dij.totalNumOfAuxVars = auxVarNumber;

% set exact optimization
pln.exactOptimization = true;