function y = matRad_interp1(xi,yi,x,extrapolation)
% interpolates 1-D data (table lookup) and utilizes griddedInterpolant
% if availabe in the used MATLAB version
%
% call
%   y = matRad_interp1(xi,yi,x)
%   y = matRad_interp1(xi,yi,x,extrapolation)
%
% input
%   xi:             sample points
%   yi:             corresponding data to sample points
%   x:              query points for interpolation
%   extrapolation:  (optional) strategy for extrapolation. Similar to
%                   interp1. NaN is the default extrapolation value
%
% output
%   y:              interpolated data
%
% Note that all input data has to be given as column vectors for a correct
% interpolation. yi can be a matrix consisting out of several 1-D datasets
% in each column, which will all be interpolated for the given query points.
%
% Reference
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent isGriddedInterpolantAvailable;

if isempty(isGriddedInterpolantAvailable)
    [env, ~] = matRad_getEnvironment();
    switch env
        case 'MATLAB'
            isGriddedInterpolantAvailable = exist('griddedInterpolant','class');
        case 'OCTAVE'
            isGriddedInterpolantAvailable = exist('griddedInterpolant');
    end
end

if nargin < 4
    extrapolation = NaN;
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
        if isscalar(extrapolation)
            y = extrapolation;
        elseif strcmp(extrapolation,'extrap')
            %In this unlikely event fall back to classic
            y = interp1(xi,yi,x,'linear',extrapolation);
        elseif strcmp(extrapolation,'nearest')
            if x < min(xi)
                y = min(yi);
            elseif x > max(xi)
                y = max(yi);
            else
                y = yi(ix1,:) + ( yi(ix2,:)-yi(ix1,:) ) * ( x - xi(ix1) ) / ( xi(ix2) - xi(ix1) );
            end
        else
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('Invalid extrapolation argument ''%s''!',extrapolation);
        end
    end

elseif isGriddedInterpolantAvailable

    if isscalar(extrapolation)
        extrapmethod = 'none';
    elseif strcmp(extrapolation,'extrap')
        extrapmethod = 'linear';
    else
        extrapmethod = extrapolation;
    end

    if size(yi,2) > 1
		% interpolation for multiple 1-D datasets
        samplePoints = {xi, 1:size(yi,2)};
        queryPoints  = {x,  1:size(yi,2)};
    else
		% interpolation for a single 1-D dataset
        samplePoints = {xi};
        queryPoints  = {x};
    end

    F = griddedInterpolant(samplePoints,yi,'linear',extrapmethod);

    y = F(queryPoints);

    if isnumeric(extrapolation) && ~isnan(extrapolation)
        y(isnan(y)) = extrapolation;
    end


else
    % for older matlab versions or octave use this code
    if isscalar(extrapolation) || strcmp(extrapolation,'extrap')
        extrapmethod = extrapolation;
    else
        matRad_cfg = MatRad_Config.instance();
        if matRad_cfg.isOctave && strcmp(extrapolation,'nearest')
            extrapmethod = NaN;
        elseif matRad_cfg.isOctave && ~strcmp(extrapolation,'linear')
            matRad_cfg.dispError('Invalid extrapolation argument ''%s''!',extrapolation);
        else
            extrapmethod = 'extrap';
        end
    end
    y = interp1(xi,yi,x(:),'linear',extrapmethod);

    %Manual nearest neighbor assignment for extrapolation
    if strcmp(extrapolation,'nearest')
        y(x(:) < min(xi),:) = repmat(min(yi),sum(x(:)<min(xi)),1);
        y(x(:) > max(xi),:) = repmat(max(yi),sum(x(:)>max(xi)),1);
    end
end
