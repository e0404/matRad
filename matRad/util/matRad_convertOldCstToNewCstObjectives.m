function newCst = matRad_convertOldCstToNewCstObjectives(cst)
% matRad function to convert cst format
% Converts a cst with struct array objectives / constraints to the new cst
% format using a cell array of objects.
% 
% call
%    newCst = matRad_convertOldCstToNewCstObjectives(cst)
%
% input
%   cst     a cst cell array that contains the old obectives as struct
%           array
%
% output 
%   newCst  copy of the input cst with all old objectives struct arrays 
%           replaced by cell arrays of Objective objects
%
% References
%   -
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Copy cst
newCst = cst;

% Loop over cst to convert objectives
for m = 1:size(cst,1)
    if ~isempty(cst{m,6})        
        %Create empty cell array in the new cst
        newCst{m,6} = cell(0);        
        %For each objective instanciate the appropriate objective object
        for n = 1:numel(cst{m,6})            
            s = matRad_DoseOptimizationFunction.convertOldOptimizationStruct(cst{m,6}(n));
            newCst{m,6}{n} = s;
        end
    end
end