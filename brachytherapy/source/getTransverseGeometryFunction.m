function GL = getTransverseGeometryFunction(obj,r)
% This function returns the value of the geometry function for the
% 1D transversal case (Rivard et al. (2004): AAPM TG-43 update, page 638f). 
% In this 1D case the recommended line source approximation is used.
%
% Output parameters:
% GL:       geometry function for the transverse plane dose distribution 
%
% Input paramaters:
% obj:      object of class Source.m
% r:        Vector containing the (radial) distance from the source center to 
%           P(r,theta), with units of cm.

% beta: is the angle subtended by the tips of the hypothetical line source 
% with respect to the calculation point, for theta_0

beta = 2*atan(obj.ActiveSourceLength./(2*r));

% Rivard et al. (2004): AAPM TG-43 update, page 638f, Eq. 4 and 11:
% GL(theta_0) = beta(theta_0)/(L*r*sin(theta_0)) with theta_0 = 90° 
GL   = beta./(obj.ActiveSourceLength*r);
        
end

