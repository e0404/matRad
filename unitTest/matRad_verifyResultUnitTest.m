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
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
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


