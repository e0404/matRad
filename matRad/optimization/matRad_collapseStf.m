function stf = matRad_collapseStf(stf,mode)
% matRad collapse stf function for simulation of 3D conformal treatments.
% Function to supress intensity-modulation for photons in order to simulate 
% 3D conformal treatments.
%
% call
%   stf = matRad_collapseStd(stf)
%
% input
%   stf:    steering information
%   mode:   collpase mode, beam or ray
%
% output
%   dijNew: collapsed dose influence matrix
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dummy collapse so that it works with sequencing

 if nargin < 2
     mode = 'beam'; % default
 end
 validatestring(mode, {'beam','ray'});

switch mode
    case 'beam'
        for i = 1:size(stf,2)
            stf(i).numOfRays = 1;
            stf(i).totalNumOfBixels = 1;
        end
    case 'ray'
        for i = 1:size(stf,2)
            stf(i).totalNumOfBixels = stf(i).numOfRays;
        end
end
