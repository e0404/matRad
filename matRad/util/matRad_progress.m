function matRad_progress(currentIndex, totalNumberOfEvaluations,scen)
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
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent isInitialized

if isempty(isInitialized)
    isInitialized = true;
end

percent = (currentIndex / totalNumberOfEvaluations) * 100;

% --- first part of the progress ---
if scen == 1
    % new line
    fprintf('\rProgress: %6.2f %%', percent);
end

% --- for scenario 2, same line ---
if scen == 2
    fprintf('   Progress: %6.2f %%', percent);
end

% --- finish line after the last angle ---
if currentIndex == totalNumberOfEvaluations && scen == 2
    fprintf('\n');
    clear isInitialized
end


