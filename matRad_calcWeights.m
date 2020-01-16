function [finalWeight, sigma_sub, posX, posY, numOfSub,X1,radius] = matRad_calcWeights(sigma_ray, n, method)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates weights for a fine sampling pencil beam
% algorithm
%
% call
%   [finalWeight, X1, sigma_sub, radius, posx, posy, numOfSub] = 
%                                  matRad_calcWeights(sigma_ray, n, method)
%
% input
%   sigma_ray:      the standard deviation of the lateral spread of the pencil
%                   beam
%   n:              number of subsample beams shells. n = 2 means we have two 
%                   lines of points around the central ray. The number of 
%                   sub-beams will be:
%                   #sb = (2*n +1)^2 for the square;
%                   #sb = (2^n -1)*6 +1 for the circle
%   method:         position shape of subsample beams. Can be 'circle' or 
%                   'square' (default 'square')
%
% output
%   finalWeight:    is the array of the weights of the sub-pencil beams. It
%                   runs over the same index as posx and posy
%   sigma_sub:      is the sigma of the gaussian of the sub-beams
%   posx & posy:    are the positions of the sub-pencil beams, returned as
%                   meshgrid-matrices if method is 'square' and as vectors 
%                   if method is 'circle'
%   numOfSub:       number of sub-pencil beams
%
% References
%   [1] http://iopscience.iop.org/0031-9155/41/8/005
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team.
%
% This file is not part of the offical matRad release
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if n~=2 && n~=3 && n~=8
    error('number of shells n not supported');
end

if ~exist('method','var')
    method = 'square';
end

if ~strcmp(method,'square') && ~strcmp(method,'circle')
    error('method not supported');
end

% This parameters come from simulations done previously
if n == 2
    if strcmp(method,'circle')
        sigma_sub = 0.8237 .* sigma_ray;
        radius    = 0.6212 .* sigma_ray;
        X1(1,:)     = 0.3866 .* sigma_ray.^2;
        X1(2,:)     = 0.6225 .* sigma_ray;
    elseif strcmp(method,'square')
        sigma_sub = 0.8409 .* sigma_ray;
        radius    = 0.5519 .* sigma_ray;
        X1(1,:)     = 0.3099 .* sigma_ray.^2;
        X1(2,:)     = 0.5556 .* sigma_ray;
    end
elseif n == 3
    if strcmp(method,'circle')
        sigma_sub = 0.7605 .* sigma_ray;
        radius    = 0.5000 .* sigma_ray;
        X1(1,:)     = 0.3006 .* sigma_ray.^2 - 1.3005 .* sigma_ray + 7.3097;
        X1(2,:)     = 0.6646 .* sigma_ray - 0.0044;
    elseif strcmp(method,'square')
        sigma_sub = 0.8409 .* sigma_ray;
        radius    = 0.5391 .* sigma_ray + 0.0856;
        X1(1,:)     = 0.3245 .* sigma_ray.^2 + 0.0001 .* sigma_ray - 0.0004;
        X1(2,:)     = 0.6290 .* sigma_ray - 0.0403;
    end
elseif n == 8 && strcmp(method,'circle')
    sigma_sub = 0.5 .* sigma_ray;
    radius    = 0.25 .* sigma_ray;
    X1(1,:)     = 0.0334 .* sigma_ray.^2 - 4.1061e-06 .* sigma_ray + 1.5047e-06;
    X1(2,:)     = 0.6 .* sigma_ray + 3.3151e-06;
end

% setting positions of sub-beams
if strcmp(method,'square')
    numOfSub = (2*n +1)^2;
    points   = linspace(-radius*n,radius*n,2*n +1);
    posX     = points'*ones(1,2*n +1);
    posY     = posX';
else
    dim = size(radius,2);
    numOfSub = (2^n -1)*6 +1;
    ang  = zeros(1,1);
    posX = zeros(1,dim);
    posY = zeros(1,dim);
    radiusShell = zeros(1,dim);
    for i = 1:n
        subsInShell = 6 * 2^(i-1);
        % this takes the sub-beams index in one shell
        ang         = cat(2, ang, pi .* linspace(0,2-2/subsInShell, subsInShell));
        radiusShell = cat(1, radiusShell, ones(subsInShell,1)*(i.*radius));
    end
    posX = cat(1, posX, bsxfun(@times,cos(ang(2:end))',radiusShell(2:end,:)));
    posY = cat(1, posY, bsxfun(@times,sin(ang(2:end))',radiusShell(2:end,:)));
end

% compute weights at positions
sig  = ones(size(posX,1),1)*X1(2,:);
normSig = ones(size(posX,1),1)*X1(1,:);

finalWeight = normSig .* (2.*pi.*sig.^2).^(-1) .* exp(-(posX.^2+posY.^2)./(2.*(sig.^2)));

