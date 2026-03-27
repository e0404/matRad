function matRad_progress(currentIndex, totalNumberOfEvaluations)
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

    percent = (currentIndex / totalNumberOfEvaluations) * 100;
    
    fprintf('\rProgress: %6.2f %%', percent);
    
    if currentIndex == totalNumberOfEvaluations
        fprintf('\n');
    end
end


