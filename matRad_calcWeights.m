function [finalWeight, X1, sigma_sub, radius, posx, posy, numOfSub]=matRad_calcWeights(sigma_ray, n, method)
% This function calculates weights for a fine sampling pencil beam
% algorithm
%
% The parameters are:
%   - sigma_ray, the standard deviation of the lateral spread of the ray
%   - n, number of subsample beams shells. n=2 means we have two lines of
%       points around the central beam. The number of sub-beams will be:
%           * #sb = (2*n +1)^2 for the square;
%           * #sb = (2^n -1)*6 +1 for the circle
%   - radius, radial distance from the central sub-beam to the nearest
%       sub-beam
%   - method, position shape of subsample beams. Can be 'circle' or 'square'
%       default 'square'
%
% Outputs are:
%   - finalWeight, is the array of the weights of the sub-gaussians. It
%       runs over the same index as posx and posy.
%   - X1, coeffincients of the weighting gaussian. In particular, if the
%      weights are a Gaussian like N*G((x,y),mu,s), so X1=[N,s];
%   - sigma_sub, is the sigma of the gaussian of the sub-beams;;
%   - posx & posy are the positions of the sub-pencil beams, returned as
%      meshgrid-matrices if method is 'square' and as vectors if method is
%      'circle'.

if n~=2 && n~=3
    error('number of shells n not supported');
end

if ~exist('method','var')
    method = 'square';
end

if ~strcmp(method,'square') && ~strcmp(method,'circle')
    error('method not supported');
end

% This parameters come from simulations done previously
if n==2
    if method=='circle'
        sigma_sub = 0.8237 * sigma_ray;
        radius = 0.6212 * sigma_ray;
        X1(1) = 0.3866 * sigma_ray^2;
        X1(2) = 0.6225 * sigma_ray;
    else % square n=2
        sigma_sub = 0.8409 * sigma_ray;
        radius = 0.5519 * sigma_ray;
        X1(1) = 0.3099 * sigma_ray^2;
        X1(2) = 0.5556 * sigma_ray;
    end
else if method=='circle'
        sigma_sub = 0.7605 * sigma_ray;
        radius = 0.5000 * sigma_ray;
        X1(1) = 0.3006 * sigma_ray^2 - 1.3005 * sigma_ray + 7.3097;
        X1(2) = 0.6646 * sigma_ray - 0.0044;
    else %square n=3
        sigma_sub = 0.8409 * sigma_ray;
        radius = 0.5391 * sigma_ray + 0.0856;
        X1(1) = 0.3245 * sigma_ray^2 + 0.0001 * sigma_ray - 0.0004;
        X1(2) = 0.6290 * sigma_ray - 0.0403;
    end
end


% setting positions of sub-beams
if strcmp(method,'square')
    numOfSub = (2*n +1)^2;
    points = linspace(-radius*(sqrt(numOfSub)-1)/2,radius*(sqrt(numOfSub)-1)/2,sqrt(numOfSub));
    posx = points'*ones(1,sqrt(numOfSub));
    posy = posx';
else
    numOfSub = (2^n -1)*6 +1;
    ang = zeros(1,1);
    posx = zeros(1,1);
    posy = zeros(1,1);
    radiusShell = zeros(1,1);
    for i=1:n
        SubsInShell = (2^i -1)*6 +1 - ((2^(i-1) -1)*6 +1 );
        % this takes the sub-beams index in one shell
        ang = cat(2, ang, pi .* linspace(0,2-2/SubsInShell, SubsInShell));
        radiusShell = cat(2, radiusShell, i.*radius.*ones(1, SubsInShell));
    end
    posx = cat(2, posx, posx(1) + radiusShell(2:end).*cos(ang(2:end)));
    posy = cat(2, posy, posy(1) + radiusShell(2:end).*sin(ang(2:end)));
end

finalWeight = zeros([1 numOfSub]);
gaussian2 = @(x, y, mux, muy ,sig) (2*pi*sig^2)^(-1) .* exp(-(x-mux).^2/(2*(sig^2)))' * exp(-(y-muy).^2/(2*(sig^2)));
for i=1:numOfSub
    finalWeight(i) = X1(1) * gaussian2(posx(i),posy(i),0,0,X1(2));
end
