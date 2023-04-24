
function F = matRad_anisotropyFunction2DInterp(r,thet,FTab)
% matRad_anisotropyFunction2DInterp anisotropy function interpolates tabulated 
%   data using interp2 ( interp technique TBD)
%   Normally called within matRad_getDoseRate(...)
%
% call
%   F = matRad_anisotropyFunction2D(r,thet,FTab)
%
% input
%   r:      array of radial distances in cm
%   thet:   array of azimuthal angles in ??
%   FTab:   tabulated consensus data of F according to the
%           following cell structure:
%           FTab{1} = AnisotropyRadialDistances
%           FTab{2} = AnisotropyPolarAngles
%           FTab{3} = AnisotropyFunctionValue
%
% output
%   F:      array of the same shape as r and thet containing the
%           interpolated and extrapolated values
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

[DataRGrid,DataThetGrid] = meshgrid(FTab{1},FTab{2});
Data(:,1) = reshape(DataRGrid,[],1);
Data(:,2) = reshape(DataThetGrid,[],1);
Value     = reshape(FTab{3},[],1);  

F = interp2(DataRGrid,DataThetGrid,FTab{3}, r, thet, 'linear');
% extrapolate for large and small values of r by taking the
% interpolation of the maximal tabulated value at this angle
% theta should be tabulated from 0?? to 180??
rmin = FTab{1}(1);
rmax = FTab{1}(end);

IndLarge = r > rmax;
IndSmall = r < rmin;
rmaxGrid = rmax*ones(sum(IndLarge),1);
rminGrid = rmin*ones(sum(IndSmall),1);
F(IndLarge) = interp2(DataRGrid,DataThetGrid,FTab{3},rmaxGrid,double(thet(IndLarge)));
F(IndSmall) = interp2(DataRGrid,DataThetGrid,FTab{3},rminGrid,double(thet(IndSmall)));

end 