function samples = matRad_getGaussianOrbitSamples(mu,SIGMA,nFrames,xr)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad orbit sampling
% 
% call
%   samples = matRad_getGaussianOrbitSamples(mu,SIGMA,nFrames,xr)
%
% input
%   mu           mean vector
%   SIGMA        covariance matrix
%   nFrames 	 number of sample frames
%   xr           if scalar, a radius for the sample. if vector, a starting coordiante
%
% output
%
% References
%   [1] http://mlss.tuebingen.mpg.de/2013/Hennig_2013_Animating_Samples_from_Gaussian_Distributions.pdf
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    samples = GPanimation(numel(mu),nFrames);
else   
    samples = GPanimation(numel(mu),nFrames,xr);
end

method = 'chol';

startMag = -10;
addMag   = 0;

switch method
    case 'chol'

        [SIGMAcol,p] = chol(SIGMA);
        if p > 0 
            disp(['Covariance Matrix not positive definite! Trying to fix.']);
            SIGMA = nearestSPD(SIGMA);
            [SIGMAcol,p] = chol(SIGMA);
            if p > 0
                error('Covariance matrix could not be fixed.');
            end
        end
        
        samples = arrayfun(@(f) mu + SIGMAcol' * samples(:,f),1:nFrames,'UniformOutput',false);
        samples = cell2mat(samples);
    case 'eig'
        SIGMAgpu = gpuArray(SIGMA);
        [V,D] = eig(SIGMAgpu);
        Q = gather(real(sqrt(complex(D)))*V);      
        samples = arrayfun(@(f) mu + Q' * samples(:,f),1:nFrames,'UniformOutput',false);
        
        samples = cell2mat(samples);
        samples = real(samples);
end
end


function X = GPanimation(d,n,xr)
% returns a matrix X of size [d,n], representing a grand circle on the
% unit d-sphere in n steps, starting at a random location. Given a kernel
% matrix K, this can be turned into a tour through the sample space, simply by
% calling chol(K)' * X; Philipp Hennig, September 2013
  
if nargin == 3
    if numel(xr) == d                      % check if we have starting coordinates
        x = xr;
        r = sqrt(sum(x.^2));
    else                                   % use radius
        r = xr;
        x = randn(d,1);
    end                              
else                                       % starting sample
    x = randn(d,1);
    r = sqrt(sum(x.^2));
end                                
 

x = x ./ sqrt(sum(x.^2));                  % project onto unit sphere
    
t = randn(d,1);                            % sample tangent direction
t = t - (t'*x) * x;                        % orthogonalise by Gram-Schmidt.
t = t ./ sqrt(sum(t.^2));                  % standardise
s = linspace(0,2*pi,n+1); s = s(1:end-1);  % space to span
t = bsxfun(@times,s,t);                    % span linspace in direction of t
X = r.* exp_map(x,t);                      % project onto sphere, re-scale 
end

function M = exp_map(mu, E)
% Computes exponential map on a sphere
D = size(E,1);
theta = sqrt(sum((E.^2)));
M = mu * cos(theta) + E .* repmat(sin(theta)./theta, D, 1);
if (any (abs (theta) <= 1e-7))
    for a = find (abs (theta) <= 1e-7)
        M (:, a) = mu;
    end % for
end % if
% M (:,abs (theta) <= 1e-7) = mu;
end % function
