function [dij,stf] = matRad_spotRemoval(dij,w,stf,thres)
% matRad postprosseing function accounting for
%       minimum number of particles per spot
%       minimum number of particles per iso-energy slice
%
% call
%   resultGUI =  matRad_postprocessing(resultGUI, dij, pln, cst, stf)

% input
%   resultGUI   struct containing optimized fluence vector
%   dij:        matRad dij struct
%   pln:        matRad pln struct
%   cst:        matRad cst struct
%   stf:        matRad stf struct
%
% output
%   resultGUI:  new w and doses in resultGUI
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2021 the matRad development team.
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

if nargin < 3 && nargout > 1
    matRad_cfg.dispError('stf is missing as input')
end

% set threshold for spot removal to 3% of the maximum weight.
if ~exist('thres','var')
    thres = 0.03;
end
newSpots = w>thres*max(w);

if ((sum(newSpots) ~= numel(w)) && sum(newSpots) ~= dij.totalNumOfBixels) && any(size(w)>1)
    dij.cutWeights = w(newSpots);
    
    dij.bixelNum = dij.bixelNum(newSpots);
    dij.rayNum = dij.rayNum(newSpots);
    dij.beamNum = dij.beamNum(newSpots);
    dij.totalNumOfBixels = sum(newSpots);
    
    dij.physicalDose{1} = dij.physicalDose{1}(:,newSpots);
    if isfield(dij,'mAlphaDose')
        dij.mAlphaDose{1} = dij.mAlphaDose{1}(:,newSpots);
        dij.mSqrtBetaDose{1} = dij.mSqrtBetaDose{1}(:,newSpots);
    end
    if isfield(dij,'mLETDose')
        dij.mLETDose{1} = dij.mLETDose{1}(:,newSpots);
    end
    [~,beamNumIdx] = unique(dij.beamNum);
    beamNumIdx = [0;beamNumIdx(2:end)-1;dij.totalNumOfBixels];
    
    for b = 1:dij.numOfBeams
        currRaysInBeam = dij.rayNum(beamNumIdx(b)+1:beamNumIdx(b+1));
        currBixelsInRay = dij.bixelNum(beamNumIdx(b)+1:beamNumIdx(b+1));
        [rayCnt,rayIdx] = unique(currRaysInBeam);
        
        if exist('stf')
            numOfBixelsPerRay = groupcounts(currRaysInBeam);
            cutRays = ismember([1:dij.numOfRaysPerBeam(b)]',rayCnt);
            if any(~cutRays)
                stf(b).ray = stf(b).ray(cutRays);
                stf(b).numOfRays = sum(cutRays);
            end
            for i = 1:stf(b).numOfRays
                bixelCurrRay{i} = currBixelsInRay(rayIdx(i):rayIdx(i)+numOfBixelsPerRay(i)-1);
            end
            for f = 1:stf(b).numOfRays
                stf(b).ray(f).energy = stf(b).ray(f).energy(bixelCurrRay{f});
                stf(b).ray(f).focusIx = stf(b).ray(f).focusIx(bixelCurrRay{f});
                stf(b).ray(f).rangeShifter = stf(b).ray(f).rangeShifter(bixelCurrRay{f});
            end
            stf(b).numOfBixelsPerRay = numOfBixelsPerRay';
            stf(b).totalNumOfBixels = sum(stf(b).numOfBixelsPerRay);
        end
        
        dij.numOfRaysPerBeam(b) = numel(rayCnt);
    end
    
    dij.totalNumOfRays = sum(dij.numOfRaysPerBeam);
    dij.numOfRemovedSpots = sum(~newSpots);
    matRad_cfg.dispWarning([num2str(sum(~newSpots)),'/',num2str(numel(newSpots)) ,' spots have been removed below ',num2str(100*thres),'%.\n'])
else
    matRad_cfg.dispWarning('no spots have been removed.')
    dij.cutWeights = w;
end
end

