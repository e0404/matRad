function [ix,radDepths] = matRad_calcRadDistsRay(ct, ...
                            isocenter, ...
                            resolution, ...
                            sourcePoint, ...
                            targetPoint, ...
                            coords)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad calculation of radiological and geometrical distances used for
% dose calcultion
% 
% call
%   [ix,radDepths] = matRad_calcRadDistsRay(ct, ...
%                       isocenter, ...
%                       resolution, ...
%                       sourcePoint, ...
%                       targetPoint, ...
%                       coords)
%
% input
%   ct:                 ct cube
%   isocenter:          isocenter
%   resolution:         resolution of the ct cube [mm]
%   sourcePoint:        source point in voxel coordinates
%   targetPoint:        target point in voxel coordinates
%   coords:             coordinates of all voxels within traced cube
%
% output
%   ix:                 indices of voxels where we want to compute dose
%                       influence data
%   radDepths:          corresponding radiological depths
%
% References
%   [1] http://www.ncbi.nlm.nih.gov/pubmed/4000088
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under 
% the terms of the Eclipse Public License 1.0 (EPL-1.0).
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.
%
% You should have received a copy of the EPL-1.0 in the file license.txt
% along with matRad. If not, see <http://opensource.org/licenses/EPL-1.0>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate radiological depths on central ray with siddon ray tracer
[alphas,l,rho,d12,ix] = matRad_siddonRayTracer(isocenter,resolution, ...
                                                sourcePoint,targetPoint, ...
                                                {ct});

% calculate geometrical distances 
geoDists = ( (coords(ix,1)-sourcePoint(1))*(targetPoint(1) - sourcePoint(1)) + ...
             (coords(ix,2)-sourcePoint(2))*(targetPoint(2) - sourcePoint(2)) + ...
             (coords(ix,3)-sourcePoint(3))*(targetPoint(3) - sourcePoint(3)) ) ...
              / d12^2;

% eq 14
% It multiply voxel intersections with \rho values.
% The zero it is neccessary for stability purpose.
d = [0 l .* rho{1}]; %Note. It is not a number "one"; it is the letter "l"

% Calculate accumulated d sum.
dCum = cumsum(d);

% This is necessary for numerical stability.
dCumIx = min([find(dCum==0,1,'last') numel(dCum)-1]);

% Calculate the radiological path
radDepths = interp1(alphas(dCumIx:end),dCum(dCumIx:end),geoDists,'linear',0);

