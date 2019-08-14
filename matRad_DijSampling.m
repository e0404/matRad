function [ixNew,bixelDoseNew] =  matRad_DijSampling(ix,bixelDose,radDepthV,rad_distancesSq,sType,Param)
% matRad dij sampling function. This function samples 
% 
% call
%   [ixNew,bixelDoseNew] =  matRad_DijSampling(ix,bixelDose,radDepthV,rad_distancesSq,r0)
%
% input
%   ix:               indices of voxels where we want to compute dose influence data
%   bixelDose:        dose at specified locations as linear vector
%   radDepthV:        radiological depth vector
%   rad_distancesSq:  squared radial distance to the central ray
%   sType:            can either be set to 'radius' or 'dose'. These are two different ways 
%                     to determine dose values that are keept as they are and dose values used for sampling
%   Param:            In the case of radius based sampling, dose values having a radial 
%                     distance below r0 [mm] are keept anyway and sampling is only done beyond r0. 
%                     In the case of dose based sampling, dose values having a relative dose greater 
%                     the threshold [0...1] are keept and sampling is done for dose values below the relative threshold  
%
% output
%   ixNew:            reduced indices of voxels where we want to compute dose influence data
%   bixelDoseNew      reduced dose at specified locations as linear vector
%
% References
%   [1] http://dx.doi.org/10.1118/1.1469633
%
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

%% define default parameters as a fallback 
defaultType                = 'radius';
deltaRadDepth              = 5;                       % step size of radiological depth
defaultLatCutOff           = 25;                      % default lateral cut off
defaultrelDoseThreshold    = 0.01;                    % default relative dose threshold

relDoseThreshold           = defaultrelDoseThreshold;
LatCutOff                  = defaultLatCutOff;
Type                       = sType;

% if the input index vector is of type logical convert it to linear indices
if islogical(ix)
   ix = find(ix); 
end

%% parse inputs
if sum(strcmp(sType,{'radius','dose'})) == 0
   Type = defaultType;
end

% if an parameter is provided then use it
if nargin>5   
    if exist('Param','var')
         if strcmp(sType,'radius')
           LatCutOff = Param;
        elseif strcmp(sType,'dose')
           relDoseThreshold = Param;
        end
    end
end

%% remember dose values inside the inner core
switch  Type
    case {'radius'}
    ixCore      = rad_distancesSq < LatCutOff^2;                 % get voxels indices having a smaller radial distance than r0
    case {'dose'}
    ixCore      = bixelDose > relDoseThreshold * max(bixelDose); % get voxels indices having a greater dose than the thresholdDose
end

bixelDoseCore       = bixelDose(ixCore);                         % save dose values that are not affected by sampling

if all(ixCore)
    %% all bixels are in the core
    %exit function with core dose only
    ixNew = ix;
    bixelDoseNew = bixelDoseCore;
else
    logIxTail           = ~ixCore;                                   % get voxels indices beyond r0
    linIxTail           = find(logIxTail);                           % convert logical index to linear index
    numTail             = numel(linIxTail);
    bixelDoseTail       = bixelDose(linIxTail);                      % dose values that are going to be reduced by sampling
    ixTail              = ix(linIxTail);                             % indices that are going to be reduced by sampling
    
    %% sample for each radiological depth the lateral halo dose
    radDepthTail        = (radDepthV(linIxTail));                    % get radiological depth in the tail
    
    % cluster radiological dephts to reduce computations
    B_r                 = int32(ceil(radDepthTail));                 % cluster radiological depths;
    maxRadDepth         = double(max(B_r));
    C                   = int32(linspace(0,maxRadDepth,round(maxRadDepth)/deltaRadDepth));     % coarse clustering of rad depths
    
    ixNew               = zeros(numTail,1);                          % inizialize new index vector
    bixelDoseNew        = zeros(numTail,1);                          % inizialize new dose vector
    linIx               = int32(1:1:numTail)';
    IxCnt               = 1;
    
    %% loop over clustered radiological depths
    for i = 1:numel(C)-1
        ixTmp              = linIx(B_r >= C(i) & B_r < C(i+1));      % extracting sub indices
        if isempty(ixTmp)
            continue
        end
        subDose            = bixelDoseTail(ixTmp);                   % get tail dose in current cluster
        subIx              = ixTail(ixTmp);                          % get indices in current cluster
        thresholdDose      = max(subDose);
        r                  = rand(numel(subDose),1);                 % get random samples
        ixSamp             = r<=(subDose/thresholdDose);
        NumSamples         = sum(ixSamp);
        
        ixNew(IxCnt:IxCnt+NumSamples-1,1)        = subIx(ixSamp);    % save new indices
        bixelDoseNew(IxCnt:IxCnt+NumSamples-1,1) = thresholdDose;    % set the dose
        IxCnt = IxCnt + NumSamples;
    end
    
    
    % cut new vectors and add inner core values
    ixNew        = [ix(ixCore);    ixNew(1:IxCnt-1)];
    bixelDoseNew = [bixelDoseCore; bixelDoseNew(1:IxCnt-1)];
end



end


