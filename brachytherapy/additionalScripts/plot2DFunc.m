% Execute from brachytherapy folder

%% plot matRad_getDoseRate2D_poly
% load brachy_HDR machine manually!

% set resolution
res = 0.05;
maxRadius = 30;

% prepare grid
x = [-maxRadius:res:maxRadius];
y = [-maxRadius:res:maxRadius];
[X,Y] = meshgrid(x,y);
[thet_rad,r] = cart2pol(X,Y);
thet = rad2deg(thet_rad);
thet(thet<0) = -thet(thet<0);
% call and plot function
DoseRate = matRad_getDoseRate2D_poly(machine,r,thet);
DoseRate(DoseRate>3) = 3;
% surf(X,Y,DoseRate,'LineStyle','none')
% xlabel('x[mm]')
% ylabel('y[mm]')
% zlabel('2D approx DoseRate')
image(DoseRate,'CDataMapping','scaled')
colorbar