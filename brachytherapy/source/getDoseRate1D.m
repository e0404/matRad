function [DoseRate] = getDoseRate1D(obj,r)
% This function returns the dose rate according to the 1D formalism 
% in Rivard et al. (2004): AAPM TG-43 update, page 639, Eq. 11.
%
% Output parameters:
% DoseRate: Dose rate in water at P(r, tetha). dose rate is returned in
%           Gy/s. 
%           When a time is specified, Dose rate is returned in Gy i.e. the
%           dose is returned !!!
%
% Input paramaters:
% obj:      object of class Source.m
% r:        Vector containing the (radial) distance from the source center to 
%           P(r,theta), with units of cm.
% time:     treatment time

global cm

if size(r,1)==1 && size(r,2)>1
    r = r(:); 
end

% TODO
% r0: The reference distance, which is 1 cm for this protocol
r0      = 1*cm;

% Sk: Air-kerma strength
Sk     = obj.SourceStrengthImplanted;
% lambda: Dose-rate constant in water (Lambda)
lambda = obj.lambda;
% GL: Geometry function G_L(r, theta_0) 
GL     = obj.getTransverseGeometryFunction(r);

% GL0: Geometry function G_L(r_0, theta_0)
GL0    = obj.getTransverseGeometryFunction(r0);

% F: 1D anisotropy function phi_an(r)
F      = obj.get1DAnisotropyFunction(r);

% gl: radial dose function g_L(r)
gl     = obj.getRadialDoseFunction(r);


% 1D dose rate formalism from Rivard et al. (2004): AAPM TG-43 update, 
% page 639, Eq. 11:
DoseRate = Sk*lambda*(GL./GL0).*gl.*F;

end

