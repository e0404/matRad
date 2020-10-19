function [cst,overlapPriorityCube] = matRad_setOverlapPriorities(cst,ctDim)
% function to handle overlap priorities 
% during fluence optimization and dose calculation. If you have overlapping 
% volumes of interest you need to inform matRad to which volume(s) the 
% intersection voxels belong.
% 
% call
%   cst = matRad_considerOverlap(cst)
%   [cst, overlapPriorityCube] = matRad_setOverlapPriorities(cst,ctDim)
%
% input
%   cst:        cst file
%   ctDim:      (optional) dimension of the ct for overlap cube claculation 
%
% output
%   cst:                updated cst file considering overlap priorities
%   overlapPriorityCube:(optional) cube visualizing the overlap priority 
%
% References
%   -
%
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

numOfCtScenarios = unique(cellfun(@(x)numel(x),cst(:,4)));

if numel(numOfCtScenarios) > 1
    error('Inconsistent number of segmentations in cst struct.');
end  

for i = 1:numOfCtScenarios
    
    % consider VOI priorities
    for j = 1:size(cst,1)
         
        idx = cst{j,4}{i};          
        
        for k = 1:size(cst,1)
            if cst{k,5}.Priority < cst{j,5}.Priority && ~(j==k) && ~isempty(cst{k,6})
                % remove indices from VOI with higher priority from current VOI
                % if an objective has been defined
                idx = setdiff(idx,cst{k,4}{i});
            end
        end
        
        cst{j,4}{i} = idx;
        
        if isempty(cst{j,4}{i}) && ~isempty(cst{j,6})
            error([cst{j,2} ': Objective(s) and/or constraints for inverse planning defined ' ...
                 'but structure overlapped by structure with higher overlap priority.' ...
                 'Objective(s) will not be considered during optimization']); 
        end
         
    end
end

%Calculate the overlap cube if requested
if nargout == 2 && nargin == 2
   overlapPriorityCube = zeros(ctDim);
    for i = 1:size(cst,1)
        overlapPriorityCube(cst{i,4}{1}) = cst{i,5}.Priority;
    end
end
    

end

