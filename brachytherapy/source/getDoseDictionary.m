function doseStruct = getDoseStruct(obj,time,distance,theta)
%GETDOSEDICTIONARY Summary of this function goes here
%   Detailed explanation goes here

if nargin == 3
    doseStruct.dictionary = obj.getDose1D(distance,time);
    doseStruct.radius     = distance;
    doseStruct.dr         = distance(2)-distance(1);
    doseStruct.theta      = [];
    doseStruct.dt         = [];

elseif nargin ==4
    [mr,mt]    = meshgrid(distance,theta);
    doseStruct.dictionary = obj.getDose2D(mr,mt);
    doseStruct.radius     = distance;
    doseStruct.dr         = distance(2)-distance(1);
    doseStruct.theta      = theta;
    doseStruct.dt         = theta(2)-theta(1);
else
    error('Not enough or too many input parameters')
end
    
    doseStruct.diameter  = obj.SourceDiameter;
    doseStruct.length    = obj.SourceLength;

end

