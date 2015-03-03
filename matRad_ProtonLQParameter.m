function [vAlpha,vBeta] = matRad_ProtonLQParameter(FreeParameter,visBool)

vAlpha = FreeParameter;
vBeta = FreeParameter;

Z = 6;
A = 12;
a = 0.0022;
p = 1.77;
E0 = (Z^2/A)^(1/p) * a^(-1/p) * R0^(1/p);

