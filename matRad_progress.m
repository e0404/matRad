function matRad_progress(CurrentIndex, TotalNumberOfEvaluations)
% This tool creates a progress bar for your loops. Call it every time your
% loop is in the next iteration and it prints your progress in % of the
% total number of steps to be carried out.
% Index has to start with '1', otherwise it doesn't work.

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
if (CurrentIndex == 1 || nargin > 2)
    fprintf('Progress: ');
end;
 
if (CurrentIndex > 1 && nargin < 3)
  Length = numel(sprintf('%3.2f %%',(CurrentIndex-1)/TotalNumberOfEvaluations*100));
  fprintf(repmat('\b',1,Length));
end
 
% Print the progress tool
fprintf('%3.2f %%',CurrentIndex/TotalNumberOfEvaluations*100);
 
% After the last iteration print a newline command
if (CurrentIndex == TotalNumberOfEvaluations)
    fprintf('\n');
end
 
end