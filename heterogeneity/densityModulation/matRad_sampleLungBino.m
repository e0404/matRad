function samples = matRad_sampleLungBino(n,p,lungDensity,numOfLungVoxels,continuous)
% matRad function for binomial sampling
%
% call
%   X = matRad_sampleBino(n,p,numOfSamples)
%
% input
%   n:              number of independent experiments
%   p:              probability (between 0 and 1)
%   numOfSamples:   number of samples for output
%
% output
%   X:              binomial samples
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% global matRad_cfg;
% matRad_cfg =  MatRad_Config.instance();

if any(~mod(n,1) == 0) && ~continuous
    error('n has to be integer for discrete distribution');
elseif ~(all(p <= 1) && all(p >= 0))
    error('p must be between 0 and 1');
elseif ~isscalar(p) && ~(numOfLungVoxels == numel(p))
    error('p array must have numOfSamples entries')
elseif ~isscalar(n) && ~(numOfLungVoxels == numel(n))
    error('n array must have numOfSamples entries')
end

% save time when testing on homogeneous phantom
if isscalar(unique(n)) && isscalar(unique(p))
    n = unique(n);
    p = unique(p);
end

if continuous
    % sample continuous beta distribution   
    a = p.*(n-1);
    b = (1-p).*(n-1);

    samples = betaincinv(rand([numOfLungVoxels,1]),a,b);
else
    % sample discrete binomial distribution
    if isscalar(n)
        samples = sum(rand([numOfLungVoxels,n]) < p, 2);
    else
        samples = zeros(numOfLungVoxels,1);
        for i = 1:numel(n)
            samples(i) = sum(rand([1,n(i)]) < p(i), 2);
        end
    end
    samples = samples ./ n;
end

samples = samples * lungDensity;

end
