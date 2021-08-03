% Execute from brachytherapy folder

%% plot matRad_getDoseRate2D_poly
% load brachy_HDR machine manually from matRad/basedata!

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