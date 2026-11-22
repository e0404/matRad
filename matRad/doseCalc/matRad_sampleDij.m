function [ixNew, bixelDoseNew] = matRad_sampleDij(ix, bixelDose, radDepthV, radDistancesSq, dijSampling, bixelWidth)
% matRad_sampleDij samples the low dose tail of a photon bixel to reduce the
% size of the dose influence matrix
%
% Voxels inside the core (defined by a lateral radius or a relative dose
% threshold) are kept. In the tail, voxels are clustered by radiological
% depth and each voxel is kept with probability dose / max(cluster dose),
% storing the cluster maximum as dose value. The expected dose is preserved.
%
% call:
%   [ixNew, bixelDoseNew] = matRad_sampleDij(ix, bixelDose, radDepthV, radDistancesSq, dijSampling)
%   [ixNew, bixelDoseNew] = matRad_sampleDij(ix, bixelDose, radDepthV, radDistancesSq, dijSampling, bixelWidth)
%
% input:
%   ix:             indices of voxels where dose influence data is computed
%   bixelDose:      dose at the specified locations as linear vector
%   radDepthV:      radiological depth vector
%   radDistancesSq: squared radial distance to the central ray
%   dijSampling:    struct with the sampling parameters relDoseThreshold,
%                   latCutOff [mm], type ('radius' or 'dose') and
%                   deltaRadDepth [mm]
%   bixelWidth:     (optional) bixel width [mm], widens the lateral core
%                   by half of the bixel diagonal
%
% output:
%   ixNew:          reduced voxel indices
%   bixelDoseNew:   reduced dose at the specified locations as linear vector
%
% References
%   [1] http://dx.doi.org/10.1118/1.1469633
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if islogical(ix)
    ix = find(ix);
end

latCutOff = dijSampling.latCutOff;
if nargin > 5 && ~isempty(bixelWidth)
    latCutOff = latCutOff + bixelWidth / sqrt(2);
end

switch dijSampling.type
    case 'radius'
        ixCore = radDistancesSq < latCutOff^2;
    case 'dose'
        ixCore = bixelDose > dijSampling.relDoseThreshold * max(bixelDose);
    otherwise
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispError('Dij Sampling mode ''%s'' not known!', dijSampling.type);
end

bixelDoseCore = bixelDose(ixCore);

if all(ixCore)
    ixNew = ix;
    bixelDoseNew = bixelDoseCore;
    return
end

linIxTail = find(~ixCore);
numTail = numel(linIxTail);
bixelDoseTail = bixelDose(linIxTail);
ixTail = ix(linIxTail);

radDepthCluster = int32(ceil(radDepthV(linIxTail)));
maxRadDepth = double(max(radDepthCluster));
numClusterEdges = max(2, floor(maxRadDepth / dijSampling.deltaRadDepth));
clusterEdges = int32(linspace(0, maxRadDepth, numClusterEdges));
clusterEdges(end) = clusterEdges(end) + 1;

ixNew = zeros(numTail, 1);
bixelDoseNew = zeros(numTail, 1);
linIx = int32(1:numTail)';
counter = 1;

for i = 1:numel(clusterEdges) - 1
    ixTmp = linIx(radDepthCluster >= clusterEdges(i) & radDepthCluster < clusterEdges(i + 1));
    if isempty(ixTmp)
        continue
    end
    subDose = bixelDoseTail(ixTmp);
    subIx = ixTail(ixTmp);
    thresholdDose = max(subDose);
    ixSamp = rand(numel(subDose), 1) <= (subDose / thresholdDose);
    numSamples = sum(ixSamp);

    ixNew(counter:counter + numSamples - 1, 1) = subIx(ixSamp);
    bixelDoseNew(counter:counter + numSamples - 1, 1) = thresholdDose;
    counter = counter + numSamples;
end

ixNew = [ix(ixCore); ixNew(1:counter - 1)];
bixelDoseNew = [bixelDoseCore; bixelDoseNew(1:counter - 1)];

end
