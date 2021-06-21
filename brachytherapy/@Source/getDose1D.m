function [Dose] = getDose1D(obj,r,time)
%GETDOSE1D Summary of this function goes here
%   Detailed explanation goes here

global cm


if size(r,1)==1 && size(r,2)>1
    r = r(:); 
end

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
Dose = Sk*lambda*(GL./GL0).*gl.*F;



if strcmp(obj.SourceType,'LDR')
    %LDR
    % decRate: decay rate
    decRate = log(2)/(obj.SourceIsotopeHalfLife);
    % DoseRate: Dose ???
    Dose = (Dose*(1-exp(-decRate*time)))./decRate;
elseif strcmp(obj.SourceType,'HDR')
    
    %    %HDR
    %    if time>0
    %        % DoseRate*time = Dose ???
    Dose = Dose*time;
else
    disp('The source must either be for HDR or LDR... LDR is assumed')
        %LDR
    % decRate: decay rate
    decRate = log(2)/(obj.SourceIsotopeHalfLife);
    % DoseRate: Dose ???
    Dose = (Dose*(1-exp(-decRate*time)))./decRate;
end

end

