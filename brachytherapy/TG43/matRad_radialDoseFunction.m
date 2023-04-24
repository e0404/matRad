function gL = matRad_radialDoseFunction(r,gLTab)
% matRad_radialDoseFuncrion interpolates tabulated data using
%   fifth order polynomial and approximates small and large distances
%   according to Rivard et al.: AAPM TG-43 update, p.669, Eq. (C1).
%   Normally called within matRad_getDoseRate(...)
%
% call
%   matRad_radialDoseFuncrion(r,gLTab)
%
% input
%   r:      array of radial distances in cm!
%   gLTab:  tabulated consensus data of gL according to the
%           following cell structure:
%           gLTab{1} = RadialDoseDistance
%           gLTab{2} = RadialDoseValue
%
% output
%   gL:     array of the same shape as r containing the interpolated
%           and extrapolated values
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

    rmin = gLTab{1}(1);
    rmax = gLTab{1}(end);
    polyCoefficients = polyfit(gLTab{1},gLTab{2},5);
    gL = zeros(size(r));
    gL(r>=rmin & r<=rmax) = polyval(polyCoefficients,r(r>=rmin & r<=rmax));
    gL(r<rmin) = gLTab{2}(1);
    gL(r>rmax) = gLTab{2}(end) + ...
                 (gLTab{2}(end)-gLTab{2}(end-1)) / (gLTab{1}(end)-...
                 gLTab{1}(end-1)).*(r(r>rmax)-gLTab{1}(end));
end
