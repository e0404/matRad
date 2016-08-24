function [ix,bixelDose] =  matRad_DijSampling(ix,bixelDose,relDoseLimits,SamplingRate)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad dij sampling function. This function samples 
% 
% call
%   [ix,bixelDose] =  matRad_DijImportanceSampling(ix,bixelDose,relDoseLimits)
%
% input
%   ix:             indices of voxels where we want to compute dose influence data
%   bixelDose:      dose at specified locations as linear vector
%   relDoseLimits:  relative dose levels used to define sampling range 
%                   e.g.[0.01 0.001] means that dij elements will be
%                   sampled in the relative dose range 1%-0.1%

%
% output
%   ix:             reduced indices of voxels where we want to compute dose influence data
%   bixelDose       reduced dose at specified locations as linear vector
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

if ~exist('SamplingRate','var')
    SamplingRate       = 0.1;
end

maxDose            = max(bixelDose);
linIxSample        = find(bixelDose < maxDose * relDoseLimits(1) & bixelDose > maxDose * relDoseLimits(2));
NumSamples         = fix(numel(linIxSample)*SamplingRate);
bixelSampDose      = bixelDose(linIxSample);
Prob               = bixelSampDose/max(bixelSampDose);
ProbNorm           = sort(Prob,'descend')./(sum(Prob));
CDF                = cumsum(ProbNorm);
randomValues       = (CDF(end)-CDF(1)).* rand(NumSamples,1) + CDF(1);
ixSamp             = interp1(CDF,linIxSample,randomValues,'nearest');
ixNew              = (bixelDose > maxDose * relDoseLimits(1)) ;
ixNew(ixSamp)      = 1;
bixelDose(ixSamp)  = maxDose * relDoseLimits(1);
bixelDose          = bixelDose(ixNew);
ix                 = ix(ixNew);




