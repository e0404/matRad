function [cst, tissueBin] = matRad_segmentationCTscan(CTdata, CTresolution, binIntervals, cst, ...
                                                      cstBodyIndex, cstTargetIndex, lungStructureName, autoSegmentLung)
% matRad segmentation of a CT scan into MCNP tissue bins
%
% 1) Classifies all voxels into the tissue bins defined by the HU intervals
%    in binIntervals (HU -> material lookup for the MCNP engine),
% 2) applies structure/morphology based tissue overrides on top of the
%    classification (see the local functions):
%      - matRad_overrideLungTissue:         lung from a contoured structure, from
%                                    auto-segmentation, or reassignment of
%                                    HU-classified lung to soft tissue
%      - matRad_overrideSkinLayer:          skin shell grown around the body hull
%      - matRad_overrideBnctTargetMaterial: boron-loaded material in a 'PTV_BNCT'
%
% NOTE: This function is private to the MCNP dose engine. It is slated to be
% replaced by a general, engine-independent HU material model with cst-based
% material overrides.
%
% call
%   [cst, tissueBin] = matRad_segmentationCTscan(CTdata, CTresolution, binIntervals, cst, ...
%                                                 cstBodyIndex, cstTargetIndex, lungStructureName, autoSegmentLung)
%
% input
%   CTdata               - CT values given in scaled HU
%   CTresolution         - CT scan resolution
%   binIntervals         - intervals for material segmentation
%   cst                  - matRad cst struct
%   cstBodyIndex         - index of the body structure in cst
%   cstTargetIndex       - index/indices of target structures in cst
%   lungStructureName    - name of a contoured lung structure (case insensitive)
%   autoSegmentLung      - segment the lung from HU if no structure is found
%
% output
%   cst                  - cst with additional LUNG/SKIN/BONE structures
%   tissueBin            - struct array with name, matIndex and voxel
%                          indices (linIndVol/matIndVol) per material,
%                          ordered as binIntervals.name
%
% Author: Lucas Sommer (Lucas.Sommer@tum.de), 06/2018
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

%% A) Classification by HU intervals given in variable binIntervals

% Important: some materials cannot be segmented using HU bin intervals
% s.th. the corresponding HUbin in the variable binIntervals is empty,
% their properties have to be defined in the last entries of binIntervals.
% For example, soft tissue with additional B-10 is defined there.
maxHUbin_nonEmpty = 1;  % Find last non-empty entry
while ~isempty(binIntervals(maxHUbin_nonEmpty + 1).HUbin) && (maxHUbin_nonEmpty <= size(binIntervals, 2))
    maxHUbin_nonEmpty = maxHUbin_nonEmpty + 1;
end

boneIndex = [];
numClassifiedVoxels = 0; % Make sure all voxels are detected
for i = 1:maxHUbin_nonEmpty
    tissueBin(i).name = binIntervals(i).name;
    tissueBin(i).matIndex = i;
    tissueBin(i).linIndVol = find((CTdata >= binIntervals(i).HUbin(1)) & CTdata < binIntervals(i).HUbin(2));

    if strcmpi(tissueBin(i).name, 'bone')
        boneIndex = i;
    end

    numClassifiedVoxels = numClassifiedVoxels + numel(tissueBin(i).linIndVol);
end

if numClassifiedVoxels ~= numel(CTdata)
    matRad_cfg.dispError('Some voxels were lost in the segmentation process!');
end

%% B) Structure/morphology based tissue overrides
airHUlimit = 300;   % Limit between air and soft tissue in scaled HU
huSmoothed = imboxfilt3(CTdata, 3, 'padding', 'replicate');   % Local mean CT values (cubic 3x3x3 box filter)

[tissueBin, cst, falseLungIdx, lungIndex] = matRad_overrideLungTissue(CTdata, tissueBin, cst, ...
                                                                      lungStructureName, autoSegmentLung, huSmoothed, airHUlimit);

[tissueBin, cst, falseLungIdx] = matRad_overrideSkinLayer(CTdata, CTresolution, tissueBin, cst, ...
                                                          cstBodyIndex, cstTargetIndex, falseLungIdx, lungIndex);

% Reassign the remaining false-lung voxels (small cavities in the body) to
% air or soft tissue according to the locally smoothed HU
if lungIndex
    for counterAdditional = 1:length(falseLungIdx)
        if huSmoothed(falseLungIdx(counterAdditional)) > airHUlimit
            tissueBin(lungIndex + 1).linIndVol(end + 1) = falseLungIdx(counterAdditional);
        else
            tissueBin(lungIndex - 1).linIndVol(end + 1) = falseLungIdx(counterAdditional);
        end
    end
    tissueBin(lungIndex - 1).linIndVol = sort(tissueBin(lungIndex - 1).linIndVol, 1);
    tissueBin(lungIndex + 1).linIndVol = sort(tissueBin(lungIndex + 1).linIndVol, 1);

    matRad_cfg.dispInfo('%s\n', ['Additional ', num2str(length(falseLungIdx)), ' voxels were redefined to be either air or soft tissue.']);
end

[tissueBin, cst] = matRad_overrideBnctTargetMaterial(tissueBin, cst);

%% C) Finalize
% Add bone to cst structure
cstCounter = size(cst, 1) + 1;
cst{cstCounter, 1} = cstCounter - 1;
cst{cstCounter, 2} = 'BONE';
cst{cstCounter, 3} = 'OAR';
cst{cstCounter, 4} = {tissueBin(boneIndex).linIndVol};
cst{cstCounter, 5}.Visible = true;
cst{cstCounter, 5}.Priority = 2;

% Make sure linear and sub-indices are the same
for counterTissueBins = 1:size(tissueBin, 2)
    tissueBin(counterTissueBins).matIndVol = [];
    if ~isempty(tissueBin(counterTissueBins).linIndVol)
        [tissueBin(counterTissueBins).matIndVol(:, 1), tissueBin(counterTissueBins).matIndVol(:, 2), ...
         tissueBin(counterTissueBins).matIndVol(:, 3)] = ind2sub(size(CTdata), tissueBin(counterTissueBins).linIndVol);
    end
end

% Check if all voxels are still assigned to a medium
numAssignedVoxels = sum(arrayfun(@(bin) length(bin.linIndVol), tissueBin));

if numAssignedVoxels > numel(CTdata)
    matRad_cfg.dispError('Something went wrong with the segmentation process. Too many voxels after segmentation.');
elseif numAssignedVoxels < numel(CTdata)
    matRad_cfg.dispError('Something went wrong with the segmentation process. Voxels were lost in segmentation process.');
else
    matRad_cfg.dispInfo('Segmentation process performed sucessfully.\n');
end

end

function [tissueBin, cst, falseLungIdx, lungIndex] = matRad_overrideLungTissue(CTdata, tissueBin, cst, ...
                                                                               lungStructureName, autoSegmentLung, huSmoothed, airHUlimit)
% Lung tissue override:
% 1) eliminate small regions falsely classified as lung by the HU intervals
%    (reassigned to air/soft tissue based on the locally smoothed HU),
% 2) replace the HU-classified lung by a contoured lung structure if
%    present, otherwise (optionally) auto-segment the lung as the largest
%    connected HU-lung region, otherwise drop the HU-lung voxels.
% Voxels removed from the lung are returned in falseLungIdx and reassigned
% by the caller after the skin layer has been cut out; lungIndex is the
% lung's tissue bin (false if the classification yielded no lung voxels).

matRad_cfg = MatRad_Config.instance();

falseLungIdx = [];

% Find material index for lung tissue from HU interval segmentation
lungIndex = find(strcmpi({tissueBin.name}, 'lung'), 1);
if isempty(lungIndex)
    matRad_cfg.dispInfo('No lung tissue from HU segmentation process.\n');
    lungIndex = false;
    return
end

if isempty(tissueBin(lungIndex).linIndVol)
    matRad_cfg.dispInfo('Segmentation via HU intervals led to zero voxels with lung tissue.\n');
    lungIndex = false;
    return
end

%% 1) Eliminate small regions that have been falsely segmented as lung
nnOfInterest = 1;   % Set size of region
minDist = sqrt(3);  % Set minimum distance to n.n. before elimination

lungMask_HUsegmentation = zeros(size(CTdata));
lungMask_HUsegmentation(tissueBin(lungIndex).linIndVol) = 1;

% Find k nearest neighbors
[a1, a2, a3] = ind2sub(size(lungMask_HUsegmentation), find(lungMask_HUsegmentation > 0));
a_ind = [a1 a2 a3];
[~, dist] = knnsearch(a_ind, a_ind, 'K', 27);

% Generate mask from tissue within lung HU interval
lungMask_processed = lungMask_HUsegmentation;

% Eliminate small pseudo-lung regions by assigning them to neighboring
% tissue intervals
eraseInd = a_ind(dist(:, nnOfInterest + 1) > minDist, :); % Number of n.n. of interest should not include voxel itself
eraseLinInd = sub2ind(size(lungMask_processed), eraseInd(:, 1), eraseInd(:, 2), eraseInd(:, 3));
lungMask_processed(eraseLinInd) = 0;

% Reassign pseudo-lung regions to tissue from surrounding HU intervals
for counterAdditional = 1:length(eraseLinInd)
    if huSmoothed(eraseLinInd(counterAdditional)) > airHUlimit
        tissueBin(lungIndex + 1).linIndVol(end + 1) = eraseLinInd(counterAdditional);
    else
        tissueBin(lungIndex - 1).linIndVol(end + 1) = eraseLinInd(counterAdditional);
    end
end

tissueBin(lungIndex - 1).linIndVol = sort(tissueBin(lungIndex - 1).linIndVol, 1);
tissueBin(lungIndex + 1).linIndVol = sort(tissueBin(lungIndex + 1).linIndVol, 1);
tissueBin(lungIndex).linIndVol = sort(find(lungMask_processed));

matRad_cfg.dispInfo('%s\n', [num2str(length(find(lungMask_HUsegmentation)) - length(find(lungMask_processed))), ...
                             ' voxels were redefined to be either air or soft tissue.']);

%% 2) Assign the lung from a contoured structure or by auto-segmentation
cstLungIndex = find(strcmpi(cst(:, 2), lungStructureName), 1);
if ~isempty(cstLungIndex)
    lungIndRegions = cst{cstLungIndex, 4}{1};
    matRad_cfg.dispInfo('Lung structure ''%s'' from the structure set used for tissue segmentation.\n', cst{cstLungIndex, 2});
elseif autoSegmentLung
    matRad_cfg.dispInfo('No lung structure ''%s'' found, lung structure will be generated...\n', lungStructureName);
    % Eliminate small regions before segmenting lung
    nnOfInterest = 9;
    eraseInd = a_ind(dist(:, nnOfInterest + 1) > minDist, :); % Number of n.n. of interest should not include voxel itself
    lungMask_processed(sub2ind(size(lungMask_processed), eraseInd(:, 1), eraseInd(:, 2), eraseInd(:, 3))) = 0;
    % Set largest connected region within body to lung
    regionConnectivity = 6; % Region connectivity (default value for 3D: 26)
    CC = bwconncomp(lungMask_processed, regionConnectivity);
    stats = regionprops3(CC, 'Volume', 'VoxelIdxList');
    [~, maxInd] = max(stats.Volume);
    lungIndRegions = stats.VoxelIdxList{maxInd, 1};
    % Generate matRad structure in cst variable
    cstCounter = size(cst, 1) + 1;
    cst{cstCounter, 1} = cstCounter - 1;
    cst{cstCounter, 2} = 'LUNG';
    cst{cstCounter, 3} = 'OAR';
    cst{cstCounter, 4} = {lungIndRegions};
    cst{cstCounter, 5}.Visible = true;
    cst{cstCounter, 5}.Priority = 2;
else
    matRad_cfg.dispInfo(['No lung structure ''%s'' found. ' ...
                         'All voxels in the HU interval for lung will be redefined as soft tissue.\n'], lungStructureName);
    matRad_cfg.dispInfo('Set the dose engine property lungStructureName or enable autoSegmentLung to segment the lung.\n');
    lungIndRegions = [];
end

% Voxels that are no longer lung are reassigned later by the caller;
% voxels that newly became lung have to leave their previous tissue bin
falseLungIdx = setdiff(tissueBin(lungIndex).linIndVol, lungIndRegions);
for binCounter = 1:size(tissueBin, 2)
    if binCounter ~= lungIndex
        tissueBin(binCounter).linIndVol = setdiff(tissueBin(binCounter).linIndVol, lungIndRegions);
    end
end
tissueBin(lungIndex).linIndVol = lungIndRegions;

end

function [tissueBin, cst, falseLungIdx] = matRad_overrideSkinLayer(CTdata, CTresolution, tissueBin, cst, ...
                                                                   cstBodyIndex, cstTargetIndex, falseLungIdx, lungIndex)
% Skin override: finds the hull of the body structure, grows a skin layer
% of skinThick mm inside it, adds the layer as an own tissue bin (and SKIN
% structure in cst) and removes its voxels from all other tissue bins as
% well as from the false-lung voxels still awaiting reassignment.

matRad_cfg = MatRad_Config.instance();

skinThick = 1; % Skin thickness [mm]

% Generate body mask
dummyMask_bodyStruct = zeros(size(CTdata));
bodyIdx = [cst{sort([cstTargetIndex, cstBodyIndex]), 4}];
bodyIdx = unique(vertcat(bodyIdx{:}));
dummyMask_bodyStruct(bodyIdx) = 1;

if skinThick > 0
    matRad_cfg.dispInfo('*****\n');
    matRad_cfg.dispInfo('%s\n', ['Skin thickness has been set to: ', num2str(skinThick), ' mm.']);
    matRad_cfg.dispInfo('*****\n');

    % Find the hull of the body structure
    linHullIdx = matRad_findBodyHull(dummyMask_bodyStruct);

    % Find k nearest neighbors of hull within body,
    % rescaled according to CT resolution s.th. origin is in voxel center
    [a1, a2, a3] = ind2sub(size(dummyMask_bodyStruct), linHullIdx);
    a_ind = [a1 * CTresolution.y - CTresolution.y / 2, ...
             a2 * CTresolution.x - CTresolution.x / 2, ...
             a3 * CTresolution.z - CTresolution.z / 2];

    [b1, b2, b3] = ind2sub(size(dummyMask_bodyStruct), find(dummyMask_bodyStruct > 0));
    b_ind = [b1 * CTresolution.y - CTresolution.y / 2, ...
             b2 * CTresolution.x - CTresolution.x / 2, ...
             b3 * CTresolution.z - CTresolution.z / 2];

    % Find neighbors and define mask
    [idx, dist] = knnsearch(b_ind, a_ind, 'K', 5000);
    dummyMask_skin = zeros(size(CTdata));

    for counterSkin = 1:size(idx, 1)
        skinInd = idx(counterSkin, dist(counterSkin, :) < skinThick);
        skinInd = b_ind(skinInd, :);
        skinInd(:, 1) = (skinInd(:, 1) + CTresolution.y / 2) / CTresolution.y;
        skinInd(:, 2) = (skinInd(:, 2) + CTresolution.x / 2) / CTresolution.x;
        skinInd(:, 3) = (skinInd(:, 3) + CTresolution.z / 2) / CTresolution.z;

        dummyMask_skin(sub2ind(size(dummyMask_bodyStruct), skinInd(:, 1), skinInd(:, 2), skinInd(:, 3))) = 1;
    end

    % Exclude presegmented air from skin segment
    airIndex = find(strcmp({tissueBin.name}, 'air'), 1);
    linSkinIdx = find(dummyMask_skin);
    idxSurfAir = intersect(linSkinIdx, tissueBin(airIndex).linIndVol);
    dummyMask_skin(idxSurfAir) = 0;

    % Add skin as tissue type
    skinBinIdx = size(tissueBin, 2) + 1;
    tissueBin(skinBinIdx).name = 'skin';
    tissueBin(skinBinIdx).matIndex = skinBinIdx;
    tissueBin(skinBinIdx).linIndVol = find(dummyMask_skin);

    % Remove skin from the false-lung voxels awaiting reassignment
    if lungIndex    % Only in case lung has been segmented using HU intervals
        falseLungIdx = setxor(falseLungIdx, intersect(falseLungIdx, tissueBin(skinBinIdx).linIndVol));
    end

    % Remove skin segment from remaining tissue segments
    for binCounter = 1:size(tissueBin, 2) - 1
        tissueBin(binCounter).linIndVol = setxor(tissueBin(binCounter).linIndVol, ...
                                                 intersect(tissueBin(binCounter).linIndVol, tissueBin(skinBinIdx).linIndVol));
    end

    % Generate matRad structure in cst variable
    cstCounter = size(cst, 1) + 1;
    cst{cstCounter, 1} = cstCounter - 1;
    cst{cstCounter, 2} = 'SKIN';
    cst{cstCounter, 3} = 'OAR';
    cst{cstCounter, 4} = {tissueBin(skinBinIdx).linIndVol};
    cst{cstCounter, 5}.Visible = true;
    cst{cstCounter, 5}.Priority = 2;
else
    % Add (empty) skin as tissue type
    skinBinIdx = size(tissueBin, 2) + 1;
    tissueBin(skinBinIdx).name = 'skin';
    tissueBin(skinBinIdx).matIndex = skinBinIdx;
    tissueBin(skinBinIdx).linIndVol = [];
    tissueBin(skinBinIdx).matIndVol = [];
end

end

function linHullIdx = matRad_findBodyHull(bodyMask)
% Finds the hull of a 3D mask by marking, slice-wise in two orientations,
% the first and last voxel of every mask row and column. Returns the linear
% indices of the hull voxels.

[yAxis, ~, zAxis] = ind2sub(size(bodyMask), find(bodyMask));
minZ = min(zAxis);
maxZ = max(zAxis);
minY = min(yAxis);
maxY = max(yAxis);

hullIdx = zeros(nnz(bodyMask), 3);
hullCounter = 1;

% Transversal slices: scan rows and columns of every connected object
for zCounter = minZ:maxZ
    stats = regionprops(bwconncomp(squeeze(bodyMask(:, :, zCounter))), 'PixelIdxList');
    for objectCounter = 1:size(stats, 1)
        dummyImage = zeros(size(squeeze(bodyMask(:, :, zCounter))));
        dummyImage(stats(objectCounter).PixelIdxList) = 1;
        for counter1 = 1:size(dummyImage, 1)
            [hullIdx, hullCounter] = matRad_markRowEnds(hullIdx, hullCounter, find(dummyImage(counter1, :)), @(i) [counter1, i, zCounter]);
        end
        for counter2 = 1:size(dummyImage, 2)
            [hullIdx, hullCounter] = matRad_markRowEnds(hullIdx, hullCounter, find(dummyImage(:, counter2)), @(i) [i, counter2, zCounter]);
        end
    end
end

% Coronal slices
for yCounter = minY:maxY
    stats = regionprops(bwconncomp(squeeze(bodyMask(yCounter, :, :))), 'PixelIdxList');
    for objectCounter = 1:size(stats, 1)
        dummyImage = zeros(size(squeeze(bodyMask(yCounter, :, :))));
        dummyImage(stats(objectCounter).PixelIdxList) = 1;
        for counter1 = 1:size(dummyImage, 1)
            [hullIdx, hullCounter] = matRad_markRowEnds(hullIdx, hullCounter, find(dummyImage(counter1, :)), @(i) [yCounter, counter1, i]);
        end
        for counter2 = 1:size(dummyImage, 2)
            [hullIdx, hullCounter] = matRad_markRowEnds(hullIdx, hullCounter, find(dummyImage(:, counter2)), @(i) [yCounter, i, counter2]);
        end
    end
end

hullIdx = hullIdx(1:hullCounter - 1, :);
linHullIdx = unique(sub2ind(size(bodyMask), hullIdx(:, 1), hullIdx(:, 2), hullIdx(:, 3)));

end

function [hullIdx, hullCounter] = matRad_markRowEnds(hullIdx, hullCounter, ind1, coordFun)
% Adds the first and last index of a mask row/column (given by coordFun) to
% the hull index list

if isempty(ind1)
    return
end

hullIdx(hullCounter, :) = coordFun(min(ind1));
hullCounter = hullCounter + 1;
if numel(ind1) > 1
    hullIdx(hullCounter, :) = coordFun(max(ind1));
    hullCounter = hullCounter + 1;
end

end

function [tissueBin, cst] = matRad_overrideBnctTargetMaterial(tissueBin, cst)
% BNCT override: if a structure named 'PTV_BNCT' is contoured, its soft
% tissue and skin voxels are moved into an own tissue bin 'bnct_material'
% (soft tissue with the B-10 loading defined in the last binIntervals
% entry); voxels of other tissues (e.g. bone) are cut from the PTV instead.

matRad_cfg = MatRad_Config.instance();

matRad_cfg.dispInfo('*****\n');
matRad_cfg.dispInfo('In case you wish to simulate BNCT irradiation make sure the PTV is called PTV_BNCT.\n');
matRad_cfg.dispInfo('*****\n');
matRad_cfg.dispInfo('Checking for PTV...\n');

bnct_cstIndex = find(strcmpi(cst(:, 2), 'PTV_BNCT'), 1);

bnctBinIdx = size(tissueBin, 2) + 1;
tissueBin(bnctBinIdx).name = 'bnct_material';
tissueBin(bnctBinIdx).matIndex = bnctBinIdx;
tissueBin(bnctBinIdx).linIndVol = [];

if isempty(bnct_cstIndex)
    matRad_cfg.dispInfo('*****\n');
    matRad_cfg.dispInfo('No PTV for BNCT detected.\n');
    matRad_cfg.dispInfo('*****\n');
    return
end

matRad_cfg.dispInfo('*****\n');
matRad_cfg.dispInfo('PTV for BNCT detected. PTV will be filled with soft tissue and B-10 density specified in segmentation variable.\n');
matRad_cfg.dispInfo('*****\n');

lindInd_PTV_BNCT = cst{bnct_cstIndex, 4}{1};
for tissueBin_counter = 1:bnctBinIdx - 1
    isSoftOrSkin = strcmpi(tissueBin(tissueBin_counter).name, 'softTissue') || ...
                   strcmpi(tissueBin(tissueBin_counter).name, 'skin');
    dummyIntersect = intersect(lindInd_PTV_BNCT, tissueBin(tissueBin_counter).linIndVol);
    if isempty(dummyIntersect)
        continue
    end

    if isSoftOrSkin
        % Move soft tissue/skin voxels from their bin into the BNCT material
        tissueBin(tissueBin_counter).linIndVol = setxor(tissueBin(tissueBin_counter).linIndVol, dummyIntersect);
        matRad_cfg.dispInfo('%s\n', [num2str(numel(dummyIntersect)), ' from ', tissueBin(tissueBin_counter).name, ...
                                     ' voxels were cut from ', tissueBin(tissueBin_counter).name, ' and associated to the BNCT material.']);
    else
        % Other tissues (e.g. bone) are cut from the PTV instead
        lindInd_PTV_BNCT = setxor(lindInd_PTV_BNCT, dummyIntersect);
        matRad_cfg.dispInfo('%s\n', [num2str(numel(dummyIntersect)), ' from ', tissueBin(tissueBin_counter).name, ' voxels were cut from PTV_BNCT.']);
    end
end

tissueBin(bnctBinIdx).linIndVol = lindInd_PTV_BNCT;

end
