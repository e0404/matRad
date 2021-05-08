function doseStruct = getDoseStruct(obj,time,distance,theta)
%GETDOSEDICTIONARY Summary of this function goes here
%   Detailed explanation goes here

distance(distance<1e-9) = 1e-9;
if nargin == 3
    if ~isempty(time)
        doseStruct.dictionary = obj.getDose1D(distance,time);
    else
        doseStruct.dictionary = obj.getDoseRate1D(distance);
    end
    doseStruct.radius     = distance;
    if length(distance)>1
    doseStruct.dr         = distance(2)-distance(1);
    end
    doseStruct.theta      = [];
    doseStruct.dt         = [];

elseif nargin ==4
   % [mr,mt]    = meshgrid(distance,theta);
    if ~isempty(time)
        doseStruct.dictionary = obj.getDose2D(mr,mt,time);
    else
        doseStruct.dictionary = abs(obj.getDoseRate2D(distance,theta));
    end
    doseStruct.radius     = distance;
    if length(distance)>1
    doseStruct.dr         = distance(2)-distance(1);
    end
    doseStruct.theta      = theta;
    if length(distance)>1
    doseStruct.dt         = theta(2)-theta(1);
    end
else
    error('Not enough or too many input parameters')
end
    
    doseStruct.sourceStrength = obj.SourceStrengthImplanted;
    doseStruct.sDiameter      = obj.SourceDiameter;
    doseStruct.sLength        = obj.SourceLength;

end

