function PhiAn = matRad_anisotropyFactor1D(r,PhiAnTab, L)
% matRad_anisotropyFactor1D anisotropy function interpolates tabulated data 
%   using fifth order polynomial and approximates small and large distances
%   according to Rivard et al.2007: Supplement to the 2004 update of the
%   AAPM Task Group No. 43 Report Eq. (2).
%   Normally called within matRad_getDoseRate(...)
%
% call
%   PhiAn = matRad_anisotropyFactor1D(r,PhiAnTab, L)
%
% input
%   r:          array of radial distances in cm!
%   PhiAnTab:   tabulated consensus data of gL according to the following
%               cell structure:
%               PhiAnTab{1} = AnisotropyFactorRadialDistance
%               PhiAnTab{2} = AnisotropyFactorValue
%
% output
%   PhiAn:      array of the same shape as r and thet containing the
%               interpolated and extrapolated values
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
    
rmin = PhiAnTab{1}(1);
rmax = PhiAnTab{1}(end);
p = polyfit(PhiAnTab{1},PhiAnTab{2},5);
PhiAn = zeros(size(r));
PhiAn(r>=rmin & r<=rmax) = polyval(p,r(r>=rmin & r<=rmax));
PhiAn(r>rmax) = PhiAnTab{2}(end);
PhiAn(r<rmin) = PhiAnTab{2}(1).*(atan(L./2./r(r<rmin))./(L.*r(r<rmin)))./(atan(L./2./rmin)./(L.*rmin));
end