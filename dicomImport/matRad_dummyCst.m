function cst = matRad_dummyCst(ct)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to create a dummy cst struct for a ct
% 
% call
%   cst = matRad_dummyCst(ct)
%
% input
%   ct: matRad ct struct
%
% output
%   cst:            matRad cst struct
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

warning('Did not find RTSS. Creating dummy segmentation for matRad.');

% allocate
cst = cell(1,6);

% fill
cst{1,1}          = 0; % first organ has number 0
cst{1,2}          = 'dummyContour';
cst{1,3}          = 'OAR';
cst{1,4}{1}       = find(ct.cube{1}>0.1);        
cst{1,5}.Priority = 1;       

% set default parameter for biological planning
cst{1,5}.alphaX   = 0.1;
cst{1,5}.betaX    = 0.05;
cst{1,5}.Visible  = 1;

% define no objcetives   
cst{1,6}          = [];

