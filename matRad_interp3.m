function y = matRad_interp3(xi,yi,zi,x,xq,yq,zq,mode,extrapVal)
% interpolates 3-D data (table lookup) 
%
% call
%   y = matRad_interp3(xi,yi,zi,x,xq,yq,zy)
%
% input
%   xi,yi,zi: grid vectors 
%   x:        data
%   xq,yq,zq: coordinates of quer points as a grid
%	
% output
%   y: interpolated data   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[env, ~] = matRad_getEnvironment();


if nargin < 8
    mode = 'linear';
end

if nargin < 9
    extrapVal = NaN;
end

switch env
    case 'MATLAB'
        y = interp3(xi,yi,zi,x,xq,yq,zq,mode,extrapVal);
    case 'OCTAVE'
        [xqMesh,yqMesh,zqMesh] = meshgrid(xq,yq,zq);
        y = interp3(xi,yi,zi,x,xqMesh,yqMesh,zqMesh,mode,extrapVal);
end

