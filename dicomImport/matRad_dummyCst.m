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

warning('Did not find RTSS. Creating dummy segmentation for matRad.');

cst = cell(1,6);

cst{1,1} = 0; % first organ has number 0
cst{1,2} = 'dummyContour';
cst{1,3} = 'OAR';
cst{1,4} = find(ct.cube>0.1);        
cst{1,5}.Priority = 1;       
% set default parameter for biological planning
cst{1,5}.TissueClass = 1; 
cst{1,5}.alphaX = 0.1;
cst{1,5}.betaX = 0.05;
cst{1,6} = []; % define no objcetives   

