function y = matRad_interp1(xi,yi,x)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolates 1-D data (table lookup) and utilizes griddedInterpolant if
% availabe in the used MATLAB version
%
% call
%   y = matRad_interp1(xi,yi,x)
%
% input
%   xi: sample points 
%   yi: corresponding data to sample points
%   x:  query points for interpolation
%	
% output
%   y: interpolated data   
%
% Note that all input data has to be given as column vectors for a correct
% interpolation. yi can be a matrix consisting out of several 1-D datasets 
% in each column, which will all be interpolated for the given query points.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent isGriddedInterpolantAvailable;

if isempty(isGriddedInterpolantAvailable)
    matRadRootDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(matRadRootDir,'tools'))
    [env, ~] = matRad_getEnvironment();
    switch env
        case 'MATLAB'
            isGriddedInterpolantAvailable = exist('griddedInterpolant','class');
        case 'OCTAVE'
            isGriddedInterpolantAvailable = exist('griddedInterpolant');
    end
end

% manual interpolation for only one query point to save time
if numel(x) == 1
    
    ix1 = find((x >= xi), 1, 'last');
    ix2 = find((x <= xi), 1, 'first');
    
    if ix1 == ix2
        y = yi(ix1,:);
    elseif ix2 == ix1 + 1
        y = yi(ix1,:) + ( yi(ix2,:)-yi(ix1,:) ) * ( x - xi(ix1) ) / ( xi(ix2) - xi(ix1) );    
    else
        error('Extrapolation not allowed');
    end
    
elseif isGriddedInterpolantAvailable
    
    if size(yi,2) > 1
		% interpolation for multiple 1-D datasets
        samplePoints = {xi, 1:size(yi,2)};
        queryPoints  = {x,  1:size(yi,2)};
    else
		% interpolation for a single 1-D dataset
        samplePoints = {xi};
        queryPoints  = {x};
    end
    
    F = griddedInterpolant(samplePoints,yi);

    y = F(queryPoints);
else
    
    % for older matlab versions use this code
    y = interp1(xi,yi,x,'linear');

end
