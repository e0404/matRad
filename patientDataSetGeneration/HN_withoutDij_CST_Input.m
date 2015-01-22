function [cst] = HN_withoutDij_CST_Input

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

%Brain Stem PRV
cst{1,3}='OAR'; cst{1,4}=1; cst{1,5}=0; cst{1,6}=1;  cst{1,7}=0;

%Brain Stem
cst{2,3}='OAR'; cst{2,4}=30; cst{2,5}=0; cst{2,6}=5;  cst{2,7}=0;

%Cerebellum
cst{3,3}='OAR'; cst{3,4}=40; cst{3,5}=0; cst{3,6}=15; cst{3,7}=0;

%Chiasma
cst{4,3}='IGNORED'; cst{4,4}=1; cst{4,5}=0; cst{4,6}=1; cst{4,7}=0;

%CTV 56
cst{5,3}='TARGET'; cst{5,4}=76; cst{5,5}=76; cst{5,6}=500; cst{5,7}=500;

%CTV 63
cst{6,3}='TARGET'; cst{6,4}=70; cst{6,5}=70; cst{6,6}=500; cst{6,7}=500;

%External
cst{7,3}='OAR'; cst{7,4}=40; cst{7,5}=0; cst{7,6}=20; cst{7,7}=0;

%GTV
cst{8,3}='OAR'; cst{8,4}=40; cst{8,5}=0; cst{8,6}=80; cst{8,7}=0;

%Larynx
cst{9,3}='OAR'; cst{9,4}=40; cst{9,5}=0; cst{9,6}=20; cst{9,7}=0;

%Lens LT
cst{10,3}='OAR'; cst{10,4}=40; cst{10,5}=0; cst{10,6}=80; cst{10,7}=0;

%Lens RT
cst{11,3}='OAR'; cst{11,4}=40; cst{11,5}=0; cst{11,6}=80; cst{11,7}=0;

%Lips
cst{12,3}='OAR'; cst{12,4}=40; cst{12,5}=0; cst{12,6}=80; cst{12,7}=0;

%Optic Nrv LT
cst{13,3}='OAR'; cst{13,4}=40; cst{13,5}=0; cst{13,6}=80; cst{13,7}=0;

%Optic Nrv RT
cst{14,3}='OAR'; cst{14,4}=40; cst{14,5}=0; cst{14,6}=80; cst{14,7}=0;

%Parotid Nrv LT
cst{15,3}='OAR'; cst{15,4}=40; cst{15,5}=0; cst{15,6}=80; cst{15,7}=0;

%Parotid Nrv RT
cst{16,3}='OAR'; cst{16,4}=40; cst{16,5}=0; cst{16,6}=80; cst{16,7}=0;

%PTV56
cst{17,3}='TARGET'; cst{17,4}=40; cst{17,5}=0; cst{17,6}=80; cst{17,7}=0;

%PTV63
cst{18,3}='TARGET'; cst{18,4}=40; cst{18,5}=0; cst{18,6}=80; cst{18,7}=0;

%PTV70
cst{19,3}='TARGET'; cst{19,4}=40; cst{19,5}=0; cst{19,6}=80; cst{19,7}=0;

%Spinal Cord
cst{20,3}='OAR'; cst{20,4}=40; cst{20,5}=0; cst{20,6}=80; cst{20,7}=0;

%Spinal Cord PRV
cst{21,3}='OAR'; cst{21,4}=40; cst{21,5}=0; cst{21,6}=80; cst{21,7}=0;

%Temp Lobe LT
cst{22,3}='OAR'; cst{22,4}=40; cst{22,5}=0; cst{22,6}=80; cst{22,7}=0;

%Temp Lobe RT
cst{23,3}='OAR'; cst{23,4}=40; cst{23,5}=0; cst{23,6}=80; cst{23,7}=0;

%TM Joint LT
cst{24,3}='IGNORED'; cst{24,4}=40; cst{24,5}=0; cst{24,6}=80; cst{24,7}=0;