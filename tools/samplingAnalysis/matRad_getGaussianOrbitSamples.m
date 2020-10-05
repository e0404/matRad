function samples = matRad_getGaussianOrbitSamples(mu,SIGMA,nFrames,varargin)
% matRad orbit sampling
% 
% call
%   samples = matRad_getGaussianOrbitSamples(mu,SIGMA,nFrames)
%   samples = matRad_getGaussianOrbitSamples(mu,SIGMA,nFrames,xr)
%   samples = matRad_getGaussianOrbitSamples(___,Name,Value)
%
% input
%   mu           mean vector
%   SIGMA        covariance matrix
%   nFrames 	 number of sample frames
%   xr           (optional) if scalar, a radius for the sample. if vector, a starting coordiante
%     
%   Optional Name-Value-Pairs:
%   Method:             can be either 'chol' (default, fast) or 'eig' 
%                       (more stable but slow)
%   MaximumTriesChol:   number of maxim tries to make matrix PSD for
%                       cholesky decomposition. Default is 10
%
%
% output
%
% References
%   [1] http://mlss.tuebingen.mpg.de/2013/Hennig_2013_Animating_Samples_from_Gaussian_Distributions.pdf
%   [2] http://www.sciencedirect.com/science/article/pii/0024379588902236
%
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

matRad_cfg = MatRad_Config.instance();


p = inputParser;

p.addRequired('mu',@(x) isvector(x) && isnumeric(x));
p.addRequired('SIGMA',@(x) size(x,1) == size(x,2) && isnumeric(x));
p.addRequired('nFrames',@(x) isscalar(x) && isnumeric(x));
p.addOptional('xr',[],@(x) isnumeric(x) && isvector(x));
p.addParameter('Method','chol',@(x) validatestring(x,{'chol'},{'eig'}));
p.addParameter('MaximumTriesChol',10,@(x) isscalar(x) && isnumeric(x));

p.parse(mu,SIGMA,nFrames,varargin{:});

method = p.Results.Method;
maxTries = p.Results.MaximumTriesChol;
xr = p.Results.xr;

%Get samples for standard Gaussian
if isempty(xr)
    samples = GPanimation(numel(mu),nFrames);
else   
    samples = GPanimation(numel(mu),nFrames,xr);
end

%Transform Samples with mean vector and covariance matrix
if strcmp(method,'chol')
    matRad_cfg.dispInfo('Creating smooth samples via cholesky decomposition...\n');
    [SIGMAcol,p] = chol(SIGMA);
    
    %If the cholesky decomposition doesn't work because the matrix is
    %not positive semi-definite, we try to fix it by finding the
    %nearest positive semidefinite matrix [2]
    nTries = 0;
    
    if p~=0
        %Equations from [2]
        [~,S,V] = svd((SIGMA + SIGMA')/2);
        SIGMA = (SIGMA+V*S*V')/2;
        SIGMA = (SIGMA + SIGMA')/2;
        
        matRad_cfg.dispInfo('Fixing Covariance matrix to be SPD: ');
        
        %To account for numerical instabilities
        while p~=0 && nTries <= maxTries
            nTries = nTries + 1;
            matRad_cfg.dispInfo('.');
            [SIGMAcol,p] = chol(SIGMA);
            if p~=0
                minEigenVal = min(eig(SIGMA));
                SIGMA = SIGMA + (-minEigenVal*nTries^2 + diag(ones(size(SIGMA,1),1)*eps(minEigenVal)));
            end
        end
        matRad_cfg.dispInfo('\n');
    end
    
    if p~=0
        matRad_cfg.dispWarning('Covariance matrix could not be fixed to be PSD in %d tries, falling back to eigenvalue method!',maxTries);
        method = 'eig';
    else 
        samples = arrayfun(@(f) mu + SIGMAcol' * samples(:,f),1:nFrames,'UniformOutput',false);
        samples = cell2mat(samples);
    end
end

if strcmp(method,'eig')
    matRad_cfg.dispInfo('Creating smooth samples via eigen decomposition...\n');
    try %Try evaluation on GPU
        SIGMAgpu = gpuArray(SIGMA);
        [V,D] = eig(SIGMAgpu);
        Q = gather(real(sqrt(complex(D)))*V);
    catch
        [V,D] = eig(SIGMA);
        Q = real(sqrt(complex(D)))*V;
    end
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
