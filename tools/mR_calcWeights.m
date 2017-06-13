function [X1, posx, posy]=mR_calcWeights(sigma_ray, n, method)

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
%   - X1, coeffincients of the weighting gaussian. In particular, if the
%      weights are a Gaussian like N*G((x,y),mu,s), so X1=[N,s];
%   - posx & posy are the positions of the sub-pencil beams, returned as
%      meshgrid-matrices if method is 'square' and as vectors if method is
%      'circle'.


startingPoint(1) = sqrt(sigma_ray^2-sigma_sub^2);
startingPoint(2) = sigma_sub^2*startingPoint(1)^2/sigma_ray^2;

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
    else
        sigma_sub = 0.8409 * sigma_ray;
        radius = 0.5519 * sigma_ray;
    end
else if method=='circle'
        sigma_sub = 0.7605 * sigma_ray;
        radius = 0.5000 * sigma_ray;
    else
        sigma_sub = 0.8409 * sigma_ray;
        radius = 0.5391 * sigma_ray + 0.0856;
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

x0 = -3*sigma_ray:sigma_ray/70:3*sigma_ray;
y0 = x0;

gaussian2 = @(x, y, mux, muy ,sig) (2*pi*sig^2)^(-1) .* exp(-(x-mux).^2/(2*(sig^2)))' * exp(-(y-muy).^2/(2*(sig^2)));
f1 =@(xi,yi) gaussian2(xi,yi,0,0,sigma_ray);
f2 = @(X,xi,yi) -f1(xi,yi);
for i=1:numOfSub
    f2 = @(X,xi,yi) f2(X,xi,yi) + X(1) .* ...
        gaussian2(posx(i),posy(i),0,0,X(2)).*gaussian2(xi,yi,posx(i),posy(i),sigma_sub);
end


f3 = @(X) 0;
x1 = -3*sigma_ray:sigma_ray/4*3:3*sigma_ray;
[xf,yf]=meshgrid(x1,x1);
numOfPoints = numel(xf);
for i=1:numOfPoints
    f3 = @(X) f3(X) + (f2(X,xf(i),yf(i)))^2;
    %             if mod(i,50)==0
    %                 disp(i)
    %             end
end

[X1,~] = fminsearch(f3, startingPoint);
todisp = [sigma_ray radius sigma_sub X1];
disp(todisp)