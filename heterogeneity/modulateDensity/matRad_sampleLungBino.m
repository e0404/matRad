function X = matRad_sampleLungBino(n,p,lungDensity,numOfLungVoxels)
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

if ~isscalar(numOfLungVoxels)
    sampleOut1 = numOfLungVoxels(1);
    sampleOut2 = numOfLungVoxels(2);
else
    sampleOut1 = numOfLungVoxels(1);
    sampleOut2 = 1;
end

if ~mod(n,1) == 0
    error('n has to be integer');
elseif ~(all(p <= 1) && all(p >= 0))
    error('p must be between 0 and 1');
elseif ~isscalar(p) && ~(sampleOut1 == numel(p))
    error('p array must have numOfSamples entries')
elseif ~isscalar(n) && ~(sampleOut1 == numel(n))
    error('n array must have numOfSamples entries')
end

if isscalar(n)
    X = sum(rand([sampleOut1,sampleOut2,n]) < p, 3);
else
    X = zeros(numel(n),1);
    for i = 1:numel(n)
        X(i) = sum(rand([sampleOut2,n(i)]) < p(i), 2);
    end
end

X = X ./ n * lungDensity;

end
