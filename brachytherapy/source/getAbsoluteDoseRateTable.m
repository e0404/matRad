function doseRateTable = getAbsoluteDoseRateTable(obj,x,y,direction)
% this function resorts the dose rates to a cartesian grid 
%   position is [0,0]
%
% Output parameters:
% doseRateTable: dose rates on a cartesian grid 
%
% Input paramaters:
% obj:          object of class Source.m
% x:            vector of coordinates in x direction
% y:            vector of coordinates in x direction
% direction:    orientation of source axis

global mm deg cm

switch nargin
    case 1
        x = [0.00	0.10	0.15	0.25	0.35	0.50	0.75	1.00	1.50	2.00	2.50	3.00	5.00	7.00]*cm;
        y = [7.00   6.00    5.00    4.00    3.00    2.50    2.00    1.50    1.00    0.75    0.50    0.25    0.10    0.00    -0.10   -0.25   -0.50   -0.75   -1.00   -1.50   -2.00   -2.50   -3.00   -4.00   -5.00   -6.00   -7.00]*cm;
        direction = [0,1];
    case 2
        y = (-7:0.1:7)*mm;
        direction = [1,0];
    case 3
        direction = [1,0];
end

% define a resoltion: 
% maxR: maximum radial distance r
maxR           = 150*mm;
% resR: resolution of radial distance r
resR           = 0.1*mm;
% resT: resolution of polar angle theta
resT           = 0.1*deg;

% radialDistance: vector of radial distances
radialDistance = 0:resR:maxR;
% theta: vector of polar angles from 0 to 90°
theta          = 0:resT:(90*deg);

% 2D dose rate
doseRate       = obj.getDoseRate2D(radialDistance,theta);

% make a cartesian grid
[mx,my] = meshgrid(x,y);

% norm: radial distance belonging to grid points (x,y) 
norm = sqrt((mx).^2+(my).^2);
% indR: indices for radial distances starting at 1. 
indR = round(norm./resR)+1;
% indT: indices for polar angles starting at 1. First calculate the angle
indT  = acos(direction(1)*mx./norm+direction(2)*my./norm);
% rearrange angles above pi/2 according to symmetry  
indT(pi    >=indT & indT>pi/2)   = pi-indT(pi>=indT & indT>pi/2);
indT(1.5*pi>=indT & indT>pi)     = indT(1.5*pi>=indT & indT>pi)-pi;
indT(2*pi  >=indT & indT>1.5*pi) = 2*pi-indT(2*pi  >=indT & indT>1.5*pi);
% get indices for angles
indT              = round(indT./resT)+1;
indT(isnan(indT)) = 1;

% Convert subscripts to linear indices
linearInd       = sub2ind(size(doseRate),indT,indR);
% resort dose rates accordingly 
doseRateTable   = doseRate(linearInd);
end

