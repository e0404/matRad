function DoseRate = getDoseRateAtPoint(obj,r,theta,time)
%GETDOSERATEATPOINT Summary of this function goes here
%   Detailed explanation goes here
tic
global mm cm deg d

switch nargin
    case 1
        r     = 10*mm;
        theta = 10*deg;
        time  = 999*d;
end

if isempty(r)
    error('The mutual distance between seed and the control point has to be defined')
end

r0       = 1*cm;
Sk       = obj.SourceStrengthImplanted;
lambda   = obj.lambda;
GL0      = obj.getTransverseGeometryFunction(r0);
gl       = obj.getRadialDoseFunction(r);


if isempty(theta)
    GL = obj.getTransverseGeometryFunction(r);
    F  = obj.get1DAnisotropyFunction(r);
else
    GL = obj.getGeometryFunction(r,theta);
    F  = obj.get2DAnisotropyFunction(r,theta);
end


DoseRate = Sk*lambda*(GL./GL0).*gl.*F;

if ~isempty(time)
    decRate = log(2)/(obj.SourceIsotopeHalfLife);
    DoseRate = (DoseRate *(1-exp(-decRate*time)))./decRate;
end
toc

end

