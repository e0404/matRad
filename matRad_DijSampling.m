function [ixNew,bixelDoseNew] =  matRad_DijSampling(ix,bixelDose,radDepthV,rad_distancesSq,r0)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad dij sampling function. This function samples 
% 
% call
%   [ixNew,bixelDoseNew] =  matRad_DijSampling(ix,bixelDose,radDepthV,rad_distancesSq,25)
%
% input
%   ix:               indices of voxels where we want to compute dose influence data
%   bixelDose:        dose at specified locations as linear vector
%   radDepthV:        radiological depth vector
%   rad_distancesSq:  squared radial distance to the central ray
%   r0:               dose values having a radial distance below r0 are keept anyway. sampling is only done beyond r0. 
%
% output
%   ixNew:            reduced indices of voxels where we want to compute dose influence data
%   bixelDoseNew      reduced dose at specified locations as linear vector
%
% References
%   [1] http://dx.doi.org/10.1118/1.1469633
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

boolFineClustering  = false;                                     % flag to switch between fine and coarse radiological depth resolution
deltaRadDepth       = 5;                                         % step size of radiological depth
LatCutOffsq         = r0^2;

%% remember dose values inside the inner core
ixCore              = rad_distancesSq < LatCutOffsq;             % get voxels indices having a smaller radial distance than r0
bixelDoseCore       = bixelDose(ixCore);                         % save dose values that are not affected by sampling

ixTail              = ~ixCore;                                   % get voxels indices beyond r0
linIxSample         = find(ixTail);                              % convert logical index to linear index
numTail             = numel(linIxSample);
bixelDoseTail       = bixelDose(linIxSample);                    % dose values that are going to be reduced by sampling
ixSampTail          = ix(linIxSample);                           % indices that are going to be reduced by sampling

%% sample for each radiological depth the lateral halo dose  
radDepthTail        = (radDepthV(linIxSample));                  % get radiological depth in the tail

% cluster radiological dephts to reduce computations
B_r                 = int32(ceil(radDepthTail));                 % cluster radiological depths; hist(B,NumOfClusters) could also be used
if boolFineClustering 
    [C,~,~]      = unique(B_r);                                                         % get unique radiological depht values == fine clustering
else 
     maxRadDepth = double(max(B_r));
     C           = int32(linspace(0,maxRadDepth,round(maxRadDepth)/deltaRadDepth));     % coarse clustering of rad depths    
end

ixNew               = zeros(numTail,1);                          % inizialize new index vector
bixelDoseNew        = zeros(numTail,1);                          % inizialize new dose vector
IxCnt               = 1;
linIx               = int32(1:1:numTail)';

% loop over clustered radiological depths
for i = 1:numel(C)-1

    ixTmp              = linIx(B_r >= C(i) & B_r < C(i+1));      % extracting sub indices
    if isempty(ixTmp)
        continue
    end 
    subDose            = bixelDoseTail(ixTmp);                   % get tail dose in current cluster
    subIx              = ixSampTail(ixTmp);                      % get indices in current cluster
    
    thresholdDose      = max(subDose); 
    r                  = rand(numel(subDose),1);                 % get random samples
    ixSamp             = r<=(subDose/thresholdDose);
    NumSamples         = sum(ixSamp);

    ixNew(IxCnt:IxCnt+NumSamples-1,1)        = subIx(ixSamp);    % save new indices
    bixelDoseNew(IxCnt:IxCnt+NumSamples-1,1) = thresholdDose;    % set the dose
    
    IxCnt = IxCnt + NumSamples;
end

% cut the vectors
ixNew        = ixNew(1:IxCnt-1);
bixelDoseNew = bixelDoseNew(1:IxCnt-1);
% add inner core values  
ixNew        = [ix(ixCore); ixNew];
bixelDoseNew = [bixelDoseCore; bixelDoseNew];

end
