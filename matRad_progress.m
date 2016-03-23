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
 
% If it's not the first step, erase the stuff printed before
if (currentIndex == 1)
    fprintf('Progress: ');
end;
 
if (currentIndex > 1)
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