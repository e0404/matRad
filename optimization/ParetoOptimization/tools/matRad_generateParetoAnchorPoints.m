function penPoints= matRad_generateParetoAnchorPoints(numObj)
% matRad function that generates the anchor points for pareto analysis
% call
%   penPoints = matRad_generateSphericalPenaltyGrid(numObj)
%
% input
%   numObj:     Number of objectives in optimization
%
% output
%   penPoints:  Matrix storing the penalty vectors
%               e.g in 3D
%                [[1,0,0]
%                 [0,1,0]
%                 [0,0,1]]
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

penPoints = eye(numObj); 

