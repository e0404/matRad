function [cst] = PROSTATE_CST_Input

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate CST column 3 to 7 for PROSTATE
% This columns are handmade generated for every organ.
% Column 3: organ class: OAR, TARGET or IGNORED. Use only Upper case.
% Column 4: maximum dose [Grays].
% Column 5: minimum dose [Grays].
% Column 6: maximum penalty.
% Column 7: minimum penalty
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

%BODY
cst{1,3}='OAR'; cst{1,4}=1; cst{1,5}=0; cst{1,6}=1;  cst{1,7}=0;

%Bladder
cst{2,3}='OAR'; cst{2,4}=30; cst{2,5}=0; cst{2,6}=5;  cst{2,7}=0;

%LT Femoral Head
cst{3,3}='OAR'; cst{3,4}=40; cst{3,5}=0; cst{3,6}=15; cst{3,7}=0;

%Lymph nodes
cst{4,3}='OAR'; cst{4,4}=1; cst{4,5}=0; cst{4,6}=1; cst{4,7}=0;

%PTV 56
cst{5,3}='TARGET'; cst{5,4}=56; cst{5,5}=56; cst{5,6}=500; cst{5,7}=500;

%PTV 68
cst{6,3}='TARGET'; cst{6,4}=68; cst{6,5}=68; cst{6,6}=500; cst{6,7}=500;

%Penile bulb
cst{7,3}='OAR'; cst{7,4}=40; cst{7,5}=0; cst{7,6}=20; cst{7,7}=0;

%Rectum
cst{8,3}='OAR'; cst{8,4}=40; cst{8,5}=0; cst{8,6}=80; cst{8,7}=0;

%Rt Femoral Head
cst{9,3}='OAR'; cst{9,4}=40; cst{9,5}=0; cst{9,6}=20; cst{9,7}=0;

%Prostate bed
cst{10,3}='TARGET'; cst{10,4}=78; cst{10,5}=78; cst{10,6}=1000; cst{10,7}=1000;