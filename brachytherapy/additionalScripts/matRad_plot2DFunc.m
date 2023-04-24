%% plot matRad_getDoseRate2D_poly
% Execute from brachytherapy folder
% load brachy_HDR machine manually from matRad/basedata!
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2021 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set resolution
res = 3;
maxRadius = 100;
DoseCutoff = 500;

% prepare grid
x = [-maxRadius:res:maxRadius];
y = [-maxRadius:res:maxRadius];
[X,Y] = meshgrid(x,y);
[thet_rad,r] = cart2pol(X,Y);
thet = rad2deg(thet_rad);
thet(thet<=0) = -thet(thet<=0);

% load basedata 
load brachy_HDR
HDRmachine = machine;
load brachy_LDR
LDRmachine = machine;
clear machine

% call and plot function
DoseRate = matRad_getDoseRate2D_poly(HDRmachine,r,thet);
DoseRate(DoseRate>DoseCutoff) = DoseCutoff;
DoseRate(DoseRate<0) = 0;

% surf(X,Y,DoseRate,'LineStyle','none')
% xlabel('x[mm]')
% ylabel('y[mm]')
% zlabel('2D approx DoseRate')
image(DoseRate,'CDataMapping','scaled')
colorbar