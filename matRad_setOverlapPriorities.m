function cst = matRad_setOverlapPriorities(cst)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to handle overlap priorities during fluence optimizaiton and
% dose calculation. If you have overlapping volumes of interest you need to
% inform matrad to which volume(s) the intersection voxels belong
% 
% call
%   cst = matRad_considerOverlap(cst)
%
% input
%   cst:    cst file
%
% output
%   cst:    updated cst file considering overlap priorities
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under 
% the terms of the Eclipse Public License 1.0 (EPL-1.0).
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.
%
% You should have received a copy of the EPL-1.0 in the file license.txt
% along with matRad. If not, see <http://opensource.org/licenses/EPL-1.0>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % consider VOI priorities
    for i = 1:size(cst,1)
         
        idx = cst{i,4};          
        
        for k = 1:size(cst,1)
            if cst{k,5}.Priority < cst{i,5}.Priority && ~(i==k)
                % remove indices from VOI with higher priority from current VOI
                idx = setdiff(idx,cst{k,4});
            end
        end
        
        cst{i,4} = idx;
        
        if isempty(cst{i,4}) && ~isempty(cst{i,6})
            warning([cst{i,2} ': Objective(s) for inverse planning defined ' ...
                 'but structure overlapped by structure with higher overlap priority.' ...
                 'Objective(s) will not be considered during optimization']); 
        end
         
    end

end

