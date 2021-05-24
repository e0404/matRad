function samples = matRad_sampleLungBino(n,p,lungDensity,numOfVoxels,continuous)
% matRad function for binomial sampling
%
% call
%   X = matRad_sampleBino(n,p,numOfSamples)
%
% input
%   n:              number of independent experiments
%   p:              probability (between 0 and 1)
%   numOfSamples:   number of samples for output, can also be an array
%                   [numOfSamples x M]
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

if nargin < 5
    continuous = false;
end

%sanity check
if ~mod(n,1) == 0
    error('n has to be integer');
elseif ~(all(p <= 1) && all(p >= 0))
    error('p must be between 0 and 1');
elseif ~isscalar(p) && ~(numOfVoxels == numel(p))
    error('p array must have numOfSamples entries')
elseif ~isscalar(n) && ~(numOfVoxels == numel(n))
    error('n array must have numOfSamples entries')
end

nVal = unique(n);
pVal = unique(p);

if continuous
    numOfDataPoints = 1e3;
    if isscalar(nVal) && isscalar(pVal)
        input = linspace(0,nVal,numOfDataPoints);
        [~,pmf,~] = matRad_continuousBino(input,nVal,pVal);
        samples = cast(discretize(rand(numOfVoxels,1), [0,pmf]),'double')/numel(input);
    else
        input = n.*linspace(0,1,numOfDataPoints);
        [~,pmf,~] = matRad_continuousBino(input,n,p);
        samples = zeros(numOfVoxels,1);
        for i = 1:numOfVoxels
            samples(i) = cast(discretize(rand(1), [0,pmf(i,:)]),'double')/numel(pmf(i,:));
        end
    end
else
    if isscalar(nVal) && isscalar(pVal)
        samples = sum(rand([numOfVoxels,1,nVal]) < pVal, 3);
    else
        samples = zeros(numel(n),1);
        for i = 1:numel(n)
            samples(i) = sum(rand([1,n(i)]) < p(i), 2);
        end
    end
    
    samples = samples ./ n;
end
samples = samples * lungDensity;
end
