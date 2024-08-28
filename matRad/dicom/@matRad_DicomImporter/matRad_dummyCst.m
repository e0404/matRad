function obj = matRad_dummyCst(obj)
% matRad function to create a dummy cst struct for a ct
% 
% call
%   obj = matRad_dummyCst(obj)
%
% input
%   ct: matRad ct struct
%
% output
%   cst: matRad cst struct
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

warning('Did not find RTSS. Creating dummy segmentation for matRad.');

% allocate
obj.cst = cell(1,6);

% fill
obj.cst{1,1}          = 0; % first organ has number 0
obj.cst{1,2}          = 'dummyContour';
obj.cst{1,3}          = 'OAR';
obj.cst{1,4}{1}       = find(obj.ct.cubeHU{1}>0.1);        
obj.cst{1,5}.Priority = 1;       

% set default parameter for biological planning
obj.cst{1,5}.alphaX   = 0.1;
obj.cst{1,5}.betaX    = 0.05;
obj.cst{1,5}.Visible  = 1;
obj.cst{1,5}.visibleColor = [0 0 0];

% define no objcetives   
obj.cst{1,6}          = [];

