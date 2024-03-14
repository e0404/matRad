function y = matRad_interp3(xi,yi,zi,x,xq,yq,zq,mode,extrapVal)
% interpolates 3-D data (table lookup) 
%
% call
%   y = matRad_interp3(xi,yi,zi,x,xq,yq,zy)
%   y = matRad_interp3(xi,yi,zi,x,xq,yq,zy,mode)
%   y = matRad_interp3(xi,yi,zi,x,xq,yq,zy,mode,extrapVal)
%
% input
%   xi,yi,zi:   grid vectors 
%   x:          data
%   xq,yq,zq:   coordinates of quer points as a grid
%   mode:       optional interpolation mode (default linear)
%   extrapVal:  (optional) value for extrapolation
%	
% output
%   y: interpolated data   
%
%   References
%     -
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
        %If we do a vector query with similar sizes only don't create a meshgrid
        if isequal(size(xq),size(yq),size(zq))
            y = interp3(xi,yi,zi,x,xq',yq',zq',mode,extrapVal);
        else
            %Here we require a meshgrid to force octave to return the correct size
            %Maybe the same thing could be achieved with a reshape?
            [xqMesh,yqMesh,zqMesh] = meshgrid(xq,yq,zq);
            y = interp3(xi,yi,zi,x,xqMesh,yqMesh,zqMesh,mode,extrapVal);
        end        
        %Old implementation
end

