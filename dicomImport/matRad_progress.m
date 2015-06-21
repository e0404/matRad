function matRad_progress(currentIndex, totalNumberOfEvaluations)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad progress bar
% 
% call
%   matRad_progress(currentIndex, totalNumberOfEvaluations)
%
% input
%   currentIndex:               current iteration index
%   totalNumberOfEvaluations:   maximum iteration index
%
% output
%   graphical display of progess. make sure there is no other output
%   written during the loop to prevent confusion
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
 
% If it's not the first step, erase the stuff printed before
if (currentIndex == 1 || nargin > 2)
    fprintf('Progress: ');
end;
 
if (currentIndex > 1 && nargin < 3)
  Length = numel(sprintf('%3.2f %%',(currentIndex-1)/totalNumberOfEvaluations*100));
  fprintf(repmat('\b',1,Length));
end
 
% Print the progress tool
fprintf('%3.2f %%',currentIndex/totalNumberOfEvaluations*100);
 
% After the last iteration print a newline command
if (currentIndex == totalNumberOfEvaluations)
    fprintf('\n');
end
 
end