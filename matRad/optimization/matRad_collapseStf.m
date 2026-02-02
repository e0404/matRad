function stf = matRad_collapseStf(stf)
% matRad collapse dij function for simulation of 3D conformal treatments.
% Function to supress intensity-modulation for photons in order to simulate 
% 3D conformal treatments.
%
% call
%   dijNew = matRad_collapseDij(dij)
%
% input
%   dij:    dose influence matrix
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

for i = 1:size(stf,2)
    stf(i).numOfRays = 1;
    stf(i).totalNumOfBixels = 1;
end