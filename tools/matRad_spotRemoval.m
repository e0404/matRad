function [dij,stf] = matRad_spotRemoval(dij,w,varargin)
% matRad spot removal tool
%
% call
%   dij =           matRad_spotRemoval(dij,w)
%   [dij,stf] =     matRad_spotRemoval(dij,w,stf)
%
% Example full call for protons
%   [dij2,stf2] = matRad_spotRemoval(dij,weights,stf)
%
% input
%   dij:                    old matRad dij struct
%   w:                      optimized matRad weights
%   varargin (optional):    stf: matRad steering struct
%                           thres: threshold for weights
%
% output
%   dij:                new matRad dij struct
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Instance of MatRad_Config class
matRad_cfg = MatRad_Config.instance();

% handle inputs
if ~isempty(varargin)
    for i = 1:nargin-2
        if isstruct(varargin{i})
            stf = varargin{i};
        elseif isscalar(varargin{i})
            thres = varargin{i};
        end
    end
end
% toggle rewriting stf if stf is an input
if exist('stf','var') && nargout > 1
    calcStf = true;
else
    calcStf = false;
end

%% calculate new spots
% For fixed threshold mode, set threshold for spot removal to 3% of the mean weight. Testing has shown differences to
% become apparent at about 5%, therefore 3% was chosen.
if ~exist('thres','var')
    thres = 0.03;
    %     thres = 0.00001;
end

% Save spots that have larger weight than the set threshold
newSpots = w>thres*mean(w);

% Generate new spots according to a percentage of the total spots
% newSpots = (1:numel(w))';
% [~,wIdx] = sort(w);
% wIdx = sort(wIdx(round(thres*numel(w)):end));
% newSpots = ismember(newSpots,wIdx);

%% rewrite dij and stf with new spots
if ((sum(newSpots) ~= numel(w)) && sum(newSpots) ~= dij.totalNumOfBixels) && any(size(w)>1)
    % save new weights
    dij.cutWeights = w(newSpots);

    % update bixel book-keeping
    dij.bixelNum = dij.bixelNum(newSpots);
    dij.rayNum = dij.rayNum(newSpots);
    dij.beamNum = dij.beamNum(newSpots);
    dij.totalNumOfBixels = sum(newSpots);

    if isfield(dij,'numParticlesPerMU') && isfield(dij,'minMU') && isfield(dij,'maxMU')
        dij.numParticlesPerMU = dij.numParticlesPerMU(newSpots);
        dij.minMU = dij.minMU(newSpots);
        dij.maxMU = dij.maxMU(newSpots);
    end

    % Freshly initialize numOfRaysPerBeam
    dij.numOfRaysPerBeam = [];

    % cut out columns in already calculated sparse matrices
    dij.physicalDose{1} = dij.physicalDose{1}(:,newSpots);
    if isfield(dij,'mAlphaDose')
        dij.mAlphaDose{1} = dij.mAlphaDose{1}(:,newSpots);
        dij.mSqrtBetaDose{1} = dij.mSqrtBetaDose{1}(:,newSpots);
    end
    if isfield(dij,'mLETDose')
        dij.mLETDose{1} = dij.mLETDose{1}(:,newSpots);
    end

    % loop through beams
    for b = 1:dij.numOfBeams

        % check if any beams have been completely removed and skip those
        if sum(dij.beamNum == b) ~= 0
            % calculate rays and indices in current beam
            currRaysInBeam = dij.rayNum(dij.beamNum == b);
            currBixelsInRay = dij.bixelNum(dij.beamNum == b);
            [rayCnt,rayIdx] = unique(currRaysInBeam);

            % save number of rays in current beam
            dij.numOfRaysPerBeam(b) = numel(rayCnt);

            % write new stf
            if calcStf
                % calculate number of bixels in each ray
                switch matRad_cfg.env
                    case 'MATLAB'
                        numOfBixelsPerRay = groupcounts(currRaysInBeam);
                    case 'OCTAVE'
                        elemts            = unique(currRaysInBeam); % MY ADDITION FOR OCTAVE
                        numOfBixelsPerRay = histc(currRaysInBeam, elemts); %MY ADDITION FOR OCTAVE
                end

                stf(b).numOfBixelsPerRay = numOfBixelsPerRay';
                stf(b).totalNumOfBixels = sum(stf(b).numOfBixelsPerRay);

                % check if any rays have been completely removed and rewrite to stf
                cutRays = ismember((1:stf(b).numOfRays)',rayCnt);
                if any(~cutRays)
                    stf(b).ray = stf(b).ray(cutRays);
                    stf(b).numOfRays = sum(cutRays);
                end

                % loop through new rays and write beam parameters for bixels that have not been removed
                for i = 1:stf(b).numOfRays
                    bixelCurrRay = currBixelsInRay(rayIdx(i):rayIdx(i)+numOfBixelsPerRay(i)-1);

                    % get fields that need to be replaced
                    fnames = fieldnames(stf(b).ray(i));
                    currFields = find(ismember(fnames,{'numParticlesPerMU','minMU','maxMU','energy','focusIx','rangeShifter'}));
                    for fieldIx = 1:numel(currFields)
                        stf(b).ray(i).(fnames{currFields(fieldIx)}) = stf(b).ray(i).(fnames{currFields(fieldIx)})(bixelCurrRay);
                    end

                end

            end
        else
            % Output warning for deleted beams and delete those beams from the new dij (and stf)
            matRad_cfg.dispWarning(['Beam ' num2str(b) ' has been deleted completely.'])

            if calcStf
                stf(b) = [];
            end

        end
    end

    % update total number of rays
    dij.totalNumOfRays = sum(dij.numOfRaysPerBeam);

    % save number of removed spots and output to console (as warning to be visible)
    dij.numOfRemovedSpots = sum(~newSpots);
    matRad_cfg.dispWarning([num2str(sum(~newSpots)),'/',num2str(numel(newSpots)) ,' spots have been removed below ',num2str(100*thres),'% of the mean weight.\n'])

    % Set MU to set minimum threshold for optimization
    dij.maxMU = Inf;
    dij.minMU = min(dij.cutWeights);
    dij.numParticlesPerMU = 1e6;
            
else
    % output warning to console
    matRad_cfg.dispWarning('no spots have been removed.')
    dij.cutWeights = w;

end
end

