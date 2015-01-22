function [cst] = TG119_CST_Input

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate CST column 3 to 7 for TG119 phantom
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
cst{1,3}='OAR'; cst{1,4}=30; cst{1,5}=0; cst{1,6}=100;  cst{1,7}=0;

%Core
cst{2,3}='OAR'; cst{2,4}=25; cst{2,5}=0; cst{2,6}=300;  cst{2,7}=0;

%Outer Target
cst{3,3}='TARGET'; cst{3,4}=50; cst{3,5}=50; cst{3,6}=1000; cst{3,7}=1000;