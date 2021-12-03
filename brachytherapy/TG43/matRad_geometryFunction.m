function GL = matRad_geometryFunction(r,thet,L)
% matRad_geometryFunction calculates 2D geometry function
%   according to Rivard et al.: AAPM TG-43, p 638 update Eq. (4)
%   Normally called within matRad_getDoseRate(...)
%
% call
%   GL = matRad_geometryFunction(r,thet,L)
%
% inputs
%   r:              array of radial distances in cm!
%   thet:           array of azimual angles in °
%   Length:         length of radiation source in cm
%
% outputs
%   GL(r,theta):    geometry function output
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2021 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% calculate solution
if thet == 90
    beta = 2*atan(L./2./r);
    GL = beta./(L.*r);
else
    GL = calcBeta(r,thet,L)./(L.*r.*sind(thet));
    GL(thet==0) = 1./(r(thet==0).^2-L^2/4);
    GL(thet==180) = 1./(r(thet==180).^2-L^2/4);
    GL(GL<0) = 0;
end

    function beta = calcBeta(r, theta,L)
        % calculate beta (see Rivard et al.: AAPM TG-43, p 637, Fig 1)
        % calculates beta from r[cm], theta [deg] and L[cm]
        % array inputs are allowed for theta

        r1 = sqrt(r.^2 + (L/2)^2 - r.*L.*cosd(180 - theta)); % cos theorem
        r2 = sqrt(r.^2 + (L/2)^2 - r.*L.*cosd(theta)); % cos theorem

        beta1 = asin(sind(180-theta).*L/2./r1); % sine theorem
        beta2 = asin(sind(theta).*L/2./r2); % sine theorem

        beta = beta1 + beta2;
    end
end