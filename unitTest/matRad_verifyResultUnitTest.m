function[Stats] = matRad_verifyResultUnitTest(Reference,Result)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad add margin function
% 
% call
%   [Stats] = matRad_verifyResultUnitTest(Reference,Result)
%
% input
%   Reference:      struct containing the reference results
%                   
%   Result:         struct containing the unit results
%
% output
%   Stats:          struct containing statistic for each evaluated field
%
% References
%   -
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

    if numel(Reference)~=numel(Result)
        warning('reference and result set are not having the same dimensions')
        Stats = 0;
        return
    end
    
    sField = fieldnames(Reference);
    
    for i=1:numel(sField)
       for j=1:length(Reference)           
           % check if field is struct - if so start recursion
           if isstruct(Reference(j).(sField{i,1}))
               Tmp(j).(sField{i,1}) = matRad_verifyResultUnitTest(Reference(j).(sField{i,1}),Result(j).(sField{i,1}));
           % check if field is vector or cube 
           elseif (isvector(Reference(j).(sField{i,1})) || numel(size((Reference(j).(sField{i,1})))) == 3 ...
                   || numel(size((Reference(j).(sField{i,1})))) == 2) && ~isstruct(Reference(j).(sField{i,1}))

                RelDiff = ((((Reference(j).(sField{i,1})+1)./(Result(j).(sField{i,1})+1))-1).*100);
                Tmp(j).(sField{i,1}).min  = min(RelDiff(:));
                Tmp(j).(sField{i,1}).max  = max(RelDiff(:));
                Tmp(j).(sField{i,1}).mean = mean(RelDiff(:));
                Tmp(j).(sField{i,1}).std  = std(RelDiff(:));
                                 
           end
       end
    end
    
    Stats = Tmp;
   
end


