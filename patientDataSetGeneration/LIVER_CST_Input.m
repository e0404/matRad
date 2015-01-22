function [cst] = LIVER_CST_Input

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate CST column 3 to 7 for LIVER
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

%CTV
cst{1,3}='TARGET'; cst{1,4}=60; cst{1,5}=60; cst{1,6}=800;  cst{1,7}=800;

%Celiac
cst{2,3}='OAR'; cst{2,4}=30; cst{2,5}=0; cst{2,6}=5;  cst{2,7}=0;

%Dose Fall Off
cst{3,3}='OAR'; cst{3,4}=40; cst{3,5}=0; cst{3,6}=15; cst{3,7}=0;

%GTV
cst{4,3}='TARGET'; cst{4,4}=76; cst{4,5}=76; cst{4,6}=1000; cst{4,7}=1000;

%Heart
cst{5,3}='OAR'; cst{5,4}=20; cst{5,5}=0; cst{5,6}=10; cst{5,7}=10;

%Kidney Left
cst{6,3}='OAR'; cst{6,4}=70; cst{6,5}=70; cst{6,6}=500; cst{6,7}=500;

%Kidney Right
cst{7,3}='OAR'; cst{7,4}=40; cst{7,5}=0; cst{7,6}=20; cst{7,7}=0;

%Large Bowel
cst{8,3}='OAR'; cst{8,4}=40; cst{8,5}=0; cst{8,6}=80; cst{8,7}=0;

%Liver
cst{9,3}='OAR'; cst{9,4}=40; cst{9,5}=0; cst{9,6}=20; cst{9,7}=0;

%PTV
cst{10,3}='TARGET'; cst{10,4}=70; cst{10,5}=70; cst{10,6}=600; cst{10,7}=600;

%SMASMV
cst{11,3}='OAR'; cst{11,4}=40; cst{11,5}=0; cst{11,6}=80; cst{11,7}=0;

%Skin
cst{12,3}='OAR'; cst{12,4}=40; cst{12,5}=0; cst{12,6}=80; cst{12,7}=0;

%Small Bowell
cst{13,3}='OAR'; cst{13,4}=40; cst{13,5}=0; cst{13,6}=80; cst{13,7}=0;

%Spinal Cord
cst{14,3}='OAR'; cst{14,4}=40; cst{14,5}=0; cst{14,6}=80; cst{14,7}=0;

%Stomach
cst{15,3}='IGNORED'; cst{15,4}=40; cst{15,5}=0; cst{15,6}=80; cst{15,7}=0;

%Duodenum
cst{16,3}='OAR'; cst{16,4}=40; cst{16,5}=0; cst{16,6}=80; cst{16,7}=0;

%Entrance
cst{17,3}='IGNORED'; cst{17,4}=40; cst{17,5}=0; cst{17,6}=80; cst{17,7}=0;

