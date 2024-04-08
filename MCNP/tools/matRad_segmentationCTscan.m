function [cst, tissueBin]=matRad_segmentationCTscan(CTdata, CTresolution, binIntervals, cst, cstBodyIndex, cstTargetIndex)
% DESCRIPTION:
% 1) Read scaled Hounsfield units and bin voxels into predefined tissue
%    bins characterized by HU intervals in binIntervals,
% 2) In case no lung has been pre-segmented define lung and reassign false
%    lung tissue to other tissue types. Furthermore, use body contour to
%    define skin layer around patient (body structure has to exist in cst).
%
% USAGE:
%       [cst, tissueBin] = segmentationCTscan(CTdata, CTresolution, binIntervals, cst, cstBodyIndex)
%
% INPUTS:
%       CTdata               - CT values given in scaled HU
%       CTresolution         - CT scan resolution
%       binIntervals         - intervals for material segementation
%       cst                  - radiotherapy structures coherent with matRad
%                              structure handling
%       cstBodyIndex         - index of body structure in cst
%
% OUTPUTS:
%       tissueBin           - structure containing all information obtained
%                             from segmentation
%       tissueBin.indices   - indices of each material specified in the bin
%                             intervals, output has the same ordering as
%                             given in binIntervals.name
%       cst                 - radiotherapy structures with additional lung
%                             and/or skin structure (depends on problem)
%
% Author: Lucas Sommer (Lucas.Sommer@tum.de), 06/2018

%% Extract information

HUbin = {binIntervals.HUbin};
tissueName = {binIntervals.name};

% Find indices in volume for each material
%% A) Segmentation by HU intervals given in variable binIntervals

dummy_indexCounter = 0; % Make sure all voxels are detected

% Important: some materials cannot be segmented using HU bin intervals
% s.th. the corresponding HUbin in the variable binIntervals is empty,
% their properties have to be defined in the last entries of binIntervals.
% For example, soft tissue with additional B-10 is defined there.
maxHUbin_nonEmpty = 1;  % Find last non-empty entry
while ~isempty(binIntervals(maxHUbin_nonEmpty+1).HUbin) && (maxHUbin_nonEmpty <= size(binIntervals,2))
    maxHUbin_nonEmpty = maxHUbin_nonEmpty +1;
end

% Segmentation using predefined HU bin intervals goes here
for i=1:maxHUbin_nonEmpty
    % Names of tissue associated to HU
    tissueBin(i).name = tissueName{i};
    
    % Tissue Indices
    tissueBin(i).matIndex = i;
    
    % Linear and matrix indices
    tissueBin(i).linIndVol = find((CTdata >= HUbin{i}(1)) & CTdata < HUbin{i}(2));   % indexing according to input
    [dum1, dum2, dum3]  = ind2sub(size(CTdata), tissueBin(i).linIndVol);
    tissueBin(i).matIndVol = [dum1, dum2, dum3];
    clear dum1; clear dum2; clear dum3;
    
    % Add bone to cst structure
    if strcmpi(tissueBin(i).name, 'bone')
        boneIndex = i;
    end
    
    dummy_indexCounter = dummy_indexCounter +numel(tissueBin(i).linIndVol);
end

if (dummy_indexCounter~=numel(CTdata))
    error('Some voxels were lost in the segmentation process!')
end

%% B) Post-processing of CT data
% 1) Eliminate small regions that have been falsely segmented as lung
% 2) Define pre-segmented lung or perform segmentation
% 3) Find body surface
% 4) Enlarge body hull for skin segmentation
% 5) Find small cavities in body and process surrounding tissue
% 6) Find BNCT PTV if present in contours

%% 1) Eliminate small regions that have been falsely segmented as lung
% Predefinitions
nnOfInterest = 1;   % Set size of region
minDist = sqrt(3);  % Set minimum distance to n.n. before elemination

airHUlimit = 300;   % Set limit for HU air value

lungIndex = 1;  % Find material index for lung tissue from HU interval segmentation
while ~strcmpi( tissueBin(lungIndex).name, 'lung') && (lungIndex<=size(tissueBin, 2))
    lungIndex = lungIndex +1;
end

if ~strcmpi( tissueBin(lungIndex).name, 'lung')  % Check if lung tissue exists in segmented material
    disp('No lung tissue from HU segmentation process.')
    lungIndex = false;
end

if ~isempty(tissueBin(lungIndex).linIndVol)
    lungMask_HUsegmentation = zeros(size(CTdata));
    lungMask_HUsegmentation(tissueBin(lungIndex).linIndVol) = 1;
else
    disp('Segmentation via HU intervals led to zero voxels with lung tissue.')
    lungIndex = false;
end


dummyCube_smoothed = imboxfilt3(CTdata,3, 'padding', 'replicate');   % Calculate local mean CT values by using cubic 3x3x3 box filter

% Only process lung material in case it has been pre-segmented
if lungIndex
    % Find k nearest neighbors
    [a1, a2, a3] = ind2sub(size(lungMask_HUsegmentation), find(lungMask_HUsegmentation>0));
    a_ind = [a1 a2 a3];
    b_ind = a_ind;
    [~, dist] = knnsearch(a_ind,b_ind,'K',27);
    
    
    % Generate mask from tissue within lung HU interval
    lungMask_processed = zeros(size(CTdata));
    lungMask_processed(tissueBin(lungIndex).linIndVol) = 1;
    
    % Eleminate small pseudo-lung regions by assigning them to neighboring
    % tissue intervals
    eraseInd = a_ind((dist(:,nnOfInterest+1)>minDist),:); % Number of n.n. of interest should not include voxel itself
    
    lungMask_processed(sub2ind(size(lungMask_processed), eraseInd(:,1), eraseInd(:,2), eraseInd(:,3))) = 0 ;
    
    % Reassign pseudo-lung regions to tissue from surrounding HU intervals
    for counterAdditional = 1:length(eraseInd)
        if dummyCube_smoothed(sub2ind(size(lungMask_processed), ...
                eraseInd(counterAdditional,1), eraseInd(counterAdditional,2), ...
                eraseInd(counterAdditional,3))) > airHUlimit
            
            tissueBin(lungIndex+1).linIndVol(end+1) = sub2ind(size(lungMask_processed), ...
                eraseInd(counterAdditional,1), ...
                eraseInd(counterAdditional,2), ...
                eraseInd(counterAdditional,3));
        elseif dummyCube_smoothed(sub2ind(size(lungMask_processed), ...
                eraseInd(counterAdditional,1), ...
                eraseInd(counterAdditional,2), eraseInd(counterAdditional,3))) <= airHUlimit
            
            tissueBin(lungIndex-1).linIndVol(end+1) = sub2ind(size(lungMask_processed), ...
                eraseInd(counterAdditional,1), ...
                eraseInd(counterAdditional,2), ...
                eraseInd(counterAdditional,3));
        end
    end
    
    tissueBin(lungIndex-1).linIndVol = sort(tissueBin(lungIndex-1).linIndVol,1);
    tissueBin(lungIndex+1).linIndVol = sort(tissueBin(lungIndex+1).linIndVol,1);
    tissueBin(lungIndex).linIndVol = sort(find(lungMask_processed));
    
    tissueBin(lungIndex-1).matIndVol = [];
    [tissueBin(lungIndex-1).matIndVol(:,1), tissueBin(lungIndex-1).matIndVol(:,2), tissueBin(lungIndex-1).matIndVol(:,3)] = ind2sub(size(CTdata), tissueBin(lungIndex-1).linIndVol);
    tissueBin(lungIndex+1).matIndVol = [];
    [tissueBin(lungIndex+1).matIndVol(:,1), tissueBin(lungIndex+1).matIndVol(:,2), tissueBin(lungIndex+1).matIndVol(:,3)] = ind2sub(size(CTdata), tissueBin(lungIndex+1).linIndVol);
    tissueBin(lungIndex).matIndVol = [];
    [tissueBin(lungIndex).matIndVol(:,1), tissueBin(lungIndex).matIndVol(:,2), tissueBin(lungIndex).matIndVol(:,3)] = ind2sub(size(CTdata), tissueBin(lungIndex).linIndVol);

    
    disp([num2str(length(find(lungMask_HUsegmentation)) - length(find(lungMask_processed))), ' voxels were redefined to be either air or soft tissue.'])
    
    % Result: lungMask_processed: Mask smoothed lung tissue from segmentation according to HU
    % intervals without small regions
    
    clear lungMask_HUsegmentation
end

%% 2. Find lung
% Overwrite lung segmentation by predefined lung structure or segment
% lung here

if lungIndex    % Check if lung voxels were generated in pre-segmentation using HU intervals
    % Check if lung is visible on CT scan
    answerLungExistance = questdlg('Is the lung visible on the CT scan?', ...
        'Lung in ROI', ...
        'Yes','No', 'No');
    
    switch answerLungExistance
        case 'Yes'  % In case lung is visible, find it
            
            % Check if lung structure is present in structure set
            cstLungIndex = 1;
            lungStructureName = {'Lung'}; % Default lung name
            
            while ~strcmpi(cst{cstLungIndex,2}, lungStructureName{1})
                cstLungIndex = cstLungIndex +1;
                if cstLungIndex > size(cst,1)
                    answer = questdlg('Is there a lung structure in your structure set (maybe it was not found)?', ...
                        'Find Pre-Segmented Lung Structure', ...
                        'Yes','No', 'No');
                    break
                end
            end
            
            
            % Handle response in case default name is not matched
            switch answer
                case 'Yes'
                    prompt = {'Please enter lung structure name:'}; % Ask for the given name
                    dlgtitle = 'Find Lung Structure';
                    lungStructureName = inputdlg(prompt,dlgtitle);
                    cstLungIndex = 1;
                    
                    while ~strcmpi(cst{cstLungIndex,2}, lungStructureName{1})   % Find lung structure
                        cstLungIndex = cstLungIndex +1;
                        if cstLungIndex > size(cst,1)
                            disp('Did not find structure. Lung structure will be generated in the following...')
                            % Eleminate small regions before segmenting lung
                            nnOfInterest = 9;
                            minDist = sqrt(3);
                            
                            eraseInd = a_ind((dist(:,nnOfInterest+1)>minDist),:); % Number of n.n. of interest should not include voxel itself
                            lungMask_processed(sub2ind(size(lungMask_processed), eraseInd(:,1), eraseInd(:,2), eraseInd(:,3))) = 0;
                            
                            % Set largest connected region within body to lung
                            regionConnectivity = 6; % Region connectivity (default value for 3D: 26)
                            CC = bwconncomp(lungMask_processed, regionConnectivity);
                            stats = regionprops3(CC,'Volume', 'VoxelIdxList');
                            
                            lungMask_final = zeros(size(CTdata));
                            [~, maxInd] = max(stats.Volume);
                            lungIndRegions = stats.VoxelIdxList{maxInd,1};
                            lungMask_final(lungIndRegions) = 1;
                            
                            falseLungIdx = setxor(tissueBin(lungIndex).linIndVol, lungIndRegions);
                            tissueBin(lungIndex).linIndVol = lungIndRegions;
                            tissueBin(lungIndex).matIndVol = [];
                            [tissueBin(lungIndex).matIndVol(:,1), tissueBin(lungIndex).matIndVol(:,2), tissueBin(lungIndex).matIndVol(:,3)] = ind2sub(size(CTdata), tissueBin(lungIndex).linIndVol);

                            break
                        else
                            falseLungIdx = setxor(tissueBin(lungIndex).linIndVol, cst{cstLungIndex,4});
                            tissueBin(lungIndex).linIndVol = cst{cstLungIndex,4};
                            tissueBin(lungIndex).matIndVol = [];
                            [tissueBin(lungIndex).matIndVol(:,1), tissueBin(lungIndex).matIndVol(:,2), tissueBin(lungIndex).matIndVol(:,3)] = ind2sub(size(CTdata), tissueBin(lungIndex).linIndVol);                            
                            
                            break
                        end
                    end
                    
                case 'No'
                    disp('Lung structure will be generated...')
                    % Eleminate small regions before segmenting lung
                    nnOfInterest = 9;
                    minDist = sqrt(3);
                    
                    eraseInd = a_ind((dist(:,nnOfInterest+1)>minDist),:); % Number of n.n. of interest should not include voxel itself
                    lungMask_processed(sub2ind(size(lungMask_processed), eraseInd(:,1), eraseInd(:,2), eraseInd(:,3))) = 0;
                    
                    % Set largest connected region within body to lung
                    regionConnectivity = 6; % Region connectivity (default value for 3D: 26)
                    CC = bwconncomp(lungMask_processed, regionConnectivity);
                    stats = regionprops3(CC,'Volume', 'VoxelIdxList');
                    
                    lungMask_final = zeros(size(CTdata));
                    [~, maxInd] = max(stats.Volume);
                    lungIndRegions = stats.VoxelIdxList{maxInd,1};
                    lungMask_final(lungIndRegions) = 1;
                    
                    falseLungIdx = setxor(tissueBin(lungIndex).linIndVol, lungIndRegions);
                    tissueBin(lungIndex).linIndVol = lungIndRegions;
                    tissueBin(lungIndex).matIndVol = [];
                    [tissueBin(lungIndex).matIndVol(:,1), tissueBin(lungIndex).matIndVol(:,2), tissueBin(lungIndex).matIndVol(:,3)] = ind2sub(size(CTdata), tissueBin(lungIndex).linIndVol);
                    
                    
                    % Generate matRad structure in cst variable
                    cstCounter = size(cst,1)+1;
                    cst{cstCounter,1} = cstCounter-1;
                    cst{cstCounter,2} = 'LUNG';
                    cst{cstCounter,3} = 'OAR';
                    cst{cstCounter,4} = {lungIndRegions};
                    cst{cstCounter,5}.Visible = true;
                    cst{cstCounter,5}.Priority = 2;

            end
            
        case 'No'   % In case there is no lung, do not do anything to find it
            falseLungIdx = tissueBin(lungIndex).linIndVol;
            tissueBin(lungIndex).linIndVol = [];
            tissueBin(lungIndex).matIndVol = [];
    end
end

% Result: lungMask_final: Mask for segmented lung

%% 3. Find body surface
% Generate body mask
dummyMask_bodyStruct = zeros(size(CTdata));
bodyIdx = [cst{sort([cstTargetIndex, cstBodyIndex]),4}];
bodyIdx = unique(vertcat(bodyIdx{:}));

dummyMask_bodyStruct(bodyIdx) = 1; % cst{cstBodyIndex,4}{1} contains body

% Find minimum and maximum z values of body axis
[yAxis,~,zAxis] = ind2sub(size(dummyMask_bodyStruct), find(dummyMask_bodyStruct));
maxZ = max(zAxis);
minZ = min(zAxis);

maxY = max(yAxis);
minY = min(yAxis);

% Find body surface
dummyMask_bodyHull = zeros(size(CTdata));

hullIdx = zeros(length(cst{cstBodyIndex,4}{1}),3);
hullCounter = 1;

for zCounter=minZ:maxZ
    stats = regionprops(bwconncomp(squeeze(dummyMask_bodyStruct(:,:,zCounter))),'Centroid', 'PixelIdxList');
    for objectCounter = 1:size(stats,1)
        dummyImage = zeros(size(squeeze(dummyMask_bodyStruct(:,:,zCounter))));
        dummyImage(stats(objectCounter).PixelIdxList)=1;
        for counter1 = 1:size(dummyImage,1)
            ind1 = find(dummyImage(counter1,:));
            if ~isempty(ind1) && numel(ind1) > 2
                hullIdx(hullCounter,:) = [counter1, min(ind1),zCounter];
                hullCounter = hullCounter +1;
                hullIdx(hullCounter,:) = [counter1, max(ind1),zCounter];
                hullCounter = hullCounter +1;
            elseif ~isempty(ind1) && numel(ind1) == 2
                hullIdx(hullCounter,:) = [counter1, min(ind1),zCounter];
                hullCounter = hullCounter +1;
                hullIdx(hullCounter,:) = [counter1, max(ind1),zCounter];
                hullCounter = hullCounter +1;
            elseif ~isempty(ind1) && numel(ind1) == 1
                hullIdx(hullCounter,:) = [counter1, min(ind1),zCounter];
            end
        end
        for counter2 = 1:size(dummyImage,2)
            ind1 = find(dummyImage(:,counter2));
            if ~isempty(ind1) && numel(ind1) > 2
                hullIdx(hullCounter,:) = [min(ind1),counter2,zCounter];
                hullCounter = hullCounter +1;
                hullIdx(hullCounter,:) = [max(ind1),counter2,zCounter];
                hullCounter = hullCounter +1;
            elseif ~isempty(ind1) && numel(ind1) == 2
                hullIdx(hullCounter,:) = [min(ind1),counter2,zCounter];
                hullCounter = hullCounter +1;
                hullIdx(hullCounter,:) = [max(ind1),counter2,zCounter];
                hullCounter = hullCounter +1;
            elseif ~isempty(ind1) && numel(ind1) == 1
                hullIdx(hullCounter,:) = [min(ind1),counter2,zCounter];
            end
            
        end
    end
end

for yCounter=minY:maxY
    stats = regionprops(bwconncomp(squeeze(dummyMask_bodyStruct(yCounter,:,:))),'Centroid', 'PixelIdxList');
    for objectCounter = 1:size(stats,1)
        dummyImage = zeros(size(squeeze(dummyMask_bodyStruct(yCounter,:,:))));
        dummyImage(stats(objectCounter).PixelIdxList)=1;
        for counter1 = 1:size(dummyImage,1)
            ind1 = find(dummyImage(counter1,:));
            if ~isempty(ind1) && numel(ind1) > 2
                hullIdx(hullCounter,:) = [yCounter, counter1, min(ind1)];
                hullCounter = hullCounter +1;
                hullIdx(hullCounter,:) = [yCounter, counter1, max(ind1)];
                hullCounter = hullCounter +1;
                
            elseif ~isempty(ind1) && numel(ind1) == 2
                hullIdx(hullCounter,:) = [yCounter, counter1, min(ind1)];
                hullCounter = hullCounter +1;
                hullIdx(hullCounter,:) = [yCounter, counter1, max(ind1)];
                hullCounter = hullCounter +1;
            elseif ~isempty(ind1) && numel(ind1) == 1
                hullIdx(hullCounter,:) = [yCounter, counter1, min(ind1)];
            end
            
        end
        for counter2 = 1:size(dummyImage,2)
            ind1 = find(dummyImage(:,counter2));
            if ~isempty(ind1) && numel(ind1) > 2
                hullIdx(hullCounter,:) = [yCounter, min(ind1), counter2];
                hullCounter = hullCounter +1;
                hullIdx(hullCounter,:) = [yCounter, max(ind1), counter2];
                hullCounter = hullCounter +1;
                
            elseif ~isempty(ind1) && numel(ind1) == 2
                hullIdx(hullCounter,:) = [yCounter, min(ind1), counter2];
                hullCounter = hullCounter +1;
                
                hullIdx(hullCounter,:) = [yCounter, max(ind1), counter2];
                hullCounter = hullCounter +1;
                
            elseif ~isempty(ind1) && numel(ind1) == 1
                hullIdx(hullCounter,:) = [yCounter, min(ind1), counter2];
                
            end
            
        end
    end
end
%
hullIdx_2 = nonzeros(hullIdx);
hullIdx_2 = reshape(hullIdx_2, [length(hullIdx_2)/3,3]);

linHullIdx = sub2ind(size(dummyMask_bodyHull), hullIdx_2(:,1), hullIdx_2(:,2), hullIdx_2(:,3));

dummyMask_bodyHull(linHullIdx) = 1;

% Result: dummyMask_bodyStruct: Mask of body structure; dummyMask_bodyHull: Mask of hull of
% body structure

%% 4. Enlarge body hull for skin segmentation

% Set skin thickness
skinThick = 1; % [mm]

if skinThick > 0
    disp('*****')
    disp(['Skin thickness has been set to: ', num2str(skinThick), ' mm.'])
    disp('*****')
    % Find k nearest neigbors of hull within body
    % Define indices of hull
    [a1, a2, a3] = ind2sub(size(dummyMask_bodyHull), find(dummyMask_bodyHull>0));
    
    a1 = a1*CTresolution.y - CTresolution.y/2;    % Rescale according to CT resolution s.th. origin is in voxel center
    a2 = a2*CTresolution.x - CTresolution.x/2;
    a3 = a3*CTresolution.z - CTresolution.z/2;
    a_ind = [a1 a2 a3];
    
    % Define body indices
    [b1, b2, b3] = ind2sub(size(dummyMask_bodyStruct), find(dummyMask_bodyStruct>0));
    
    b1 = b1*CTresolution.y - CTresolution.y/2;    % Rescale according to CT resolution
    b2 = b2*CTresolution.x - CTresolution.x/2;
    b3 = b3*CTresolution.z - CTresolution.z/2;
    b_ind = [b1 b2 b3];
    
    % Find neighbors and define mask
    [idx, dist] = knnsearch(b_ind, a_ind,'K',5000);
    dummyMask_skin = zeros(size(CTdata));
    
    for counterSkin = 1:size(idx,1)
        skinInd = idx(counterSkin,dist(counterSkin,:)<skinThick);
        skinInd = b_ind(skinInd, :);
        skinInd(:,1) = (skinInd(:,1) + CTresolution.y/2)/CTresolution.y;
        skinInd(:,2) = (skinInd(:,2) + CTresolution.x/2)/CTresolution.x;
        skinInd(:,3) = (skinInd(:,3) + CTresolution.z/2)/CTresolution.z;
        
        dummyMask_skin(sub2ind(size(dummyMask_bodyStruct), skinInd(:,1), skinInd(:,2), skinInd(:,3))) = 1;
        clear skinInd;
    end
    
    % Exclude presegmented air from skin segment
    airIndex = 1;  % Find material index for lung tissue from HU interval segmentation
    while ~strcmp( tissueBin(airIndex).name, 'air')
        airIndex = airIndex +1;
    end
    
    linSkinIdx = find(dummyMask_skin);
    idxSurfAir = intersect(linSkinIdx, tissueBin(airIndex).linIndVol);
    dummyMask_skin(idxSurfAir) = 0;
    
    % Add skin as tissue type
    skinBinIdx = size(tissueBin,2)+1;
    tissueBin(skinBinIdx).name = 'skin';
    tissueBin(skinBinIdx).matIndex = skinBinIdx;
    tissueBin(skinBinIdx).linIndVol = find(dummyMask_skin);
    
    [dumIdx1, dumIdx2, dumIdx3] = ind2sub(size(dummyMask_skin), tissueBin(skinBinIdx).linIndVol);
    tissueBin(skinBinIdx).matIndVol = [dumIdx1 dumIdx2 dumIdx3];
    
    % Remove skin from other tissue bins
    % Start with voxel uncorrectly assigend to lung
    if lungIndex    % Only in case lung has been segmented using HU intervals
        falseLungIdx = setxor(falseLungIdx, intersect(falseLungIdx, tissueBin(skinBinIdx).linIndVol));
    end
    
    % Remove skin segment from remaining tissue segments
    for binCounter = 1:size(tissueBin,2)-1
        tissueBin(binCounter).linIndVol = setxor(tissueBin(binCounter).linIndVol, intersect(tissueBin(binCounter).linIndVol, tissueBin(skinBinIdx).linIndVol));
        [dumIdx1, dumIdx2, dumIdx3] = ind2sub(size(dummyMask_skin), tissueBin(binCounter).linIndVol);
        tissueBin(binCounter).matIndVol = [dumIdx1, dumIdx2, dumIdx3];
    end
    
    % Generate matRad structure in cst variable
    cstCounter = size(cst,1)+1;
    cst{cstCounter,1} = cstCounter-1;
    cst{cstCounter,2} = 'SKIN';
    cst{cstCounter,3} = 'OAR';
    cst{cstCounter,4} = {tissueBin(skinBinIdx).linIndVol};
    cst{cstCounter,5}.Visible = true;
    cst{cstCounter,5}.Priority = 2;
    
elseif skinThick == 0
    % Add skin as tissue type
    skinBinIdx = size(tissueBin,2)+1;
    tissueBin(skinBinIdx).name = 'skin';
    tissueBin(skinBinIdx).matIndex = skinBinIdx;
    tissueBin(skinBinIdx).linIndVol = [];
    tissueBin(skinBinIdx).matIndVol = [];
end

% Result: dummyMask_skin: Mask of skin structure


%% 5. Find small cavities in body and process surrounding tissue

if lungIndex    % Only in case lung has been segmented using HU intervals
    falseLungIdx = intersect(falseLungIdx, setxor(falseLungIdx, tissueBin(skinBinIdx).linIndVol));
    dummyMask_falseLungToProcess = zeros(size(CTdata));
    dummyMask_falseLungToProcess(falseLungIdx) = 1;
    
    % Reassign pseudo-lung regions to another tissue
    for counterAdditional = 1:length(falseLungIdx)
        if dummyCube_smoothed(falseLungIdx(counterAdditional)) > airHUlimit
            tissueBin(lungIndex+1).linIndVol(end+1) = falseLungIdx(counterAdditional);                      
        elseif dummyCube_smoothed(falseLungIdx(counterAdditional)) <= airHUlimit
            tissueBin(lungIndex-1).linIndVol(end+1) = falseLungIdx(counterAdditional);
        end
    end
    
    tissueBin(lungIndex-1).linIndVol = sort(tissueBin(lungIndex-1).linIndVol,1);
    tissueBin(lungIndex-1).matIndVol = [];
    [tissueBin(lungIndex-1).matIndVol(:,1), tissueBin(lungIndex-1).matIndVol(:,2), tissueBin(lungIndex-1).matIndVol(:,3)] = ind2sub(size(CTdata), tissueBin(lungIndex-1).linIndVol);
    
    
    tissueBin(lungIndex+1).linIndVol = sort(tissueBin(lungIndex+1).linIndVol,1);
    tissueBin(lungIndex+1).matIndVol = [];
    [tissueBin(lungIndex+1).matIndVol(:,1), tissueBin(lungIndex+1).matIndVol(:,2), tissueBin(lungIndex+1).matIndVol(:,3)] = ind2sub(size(CTdata), tissueBin(lungIndex+1).linIndVol);
    
    
    disp(['Additional ', num2str(length(falseLungIdx)), ' voxels were redefined to be either air or soft tissue.'])
    
end

%% 6. Find BNCT 
disp('*****')
disp('In case you wish to simulate BNCT irradiation make sure the PTV is called PTV_BNCT.')
disp('*****')
disp('Checking for PTV...')


% Find BNCT PTV
for cstCounter = 1:size(cst,1)
    if strcmpi(cst{cstCounter,2}, 'PTV_BNCT')
        disp('*****')
        disp('PTV for BNCT detected. PTV will be filled with soft tissu and B-10 density specified in segmentation variable.')
        disp('*****')
        
        bnct_Control = true;
        bnct_cstIndex = cstCounter;
    else
        continue
    end
end

% Process BNCT volume and generate specific material with B-10 for MCNP
% simulation
if exist('bnct_Control') && bnct_Control
    lindInd_PTV_BNCT = cst{bnct_cstIndex,4};
    for tissueBin_counter = 1:size(tissueBin,2)
        if ~strcmpi(tissueBin(tissueBin_counter).name, 'softTissue') && ~strcmpi(tissueBin(tissueBin_counter).name, 'skin')
            dummyIntersect = intersect(lindInd_PTV_BNCT{1}, tissueBin(tissueBin_counter).linIndVol);
            if ~isempty(dummyIntersect)
                lindInd_PTV_BNCT{1} = setxor(lindInd_PTV_BNCT{1}, dummyIntersect);
                disp('*****')
                disp([num2str(numel(dummyIntersect)), ' from ', tissueBin(tissueBin_counter).name,' voxels were cut from PTV_BNCT.'])
                disp('*****')
            end
        elseif strcmpi(tissueBin(tissueBin_counter).name, 'softTissue') || strcmpi(tissueBin(tissueBin_counter).name, 'skin')
            dummyIntersect = intersect(lindInd_PTV_BNCT{1}, tissueBin(tissueBin_counter).linIndVol);
            if ~isempty(dummyIntersect)
                tissueBin(tissueBin_counter).linIndVol = setxor(tissueBin(tissueBin_counter).linIndVol, dummyIntersect);
                disp('*****')
                disp([num2str(numel(dummyIntersect)), ' from ', tissueBin(tissueBin_counter).name, ' voxels were cut from ', tissueBin(tissueBin_counter).name, ' and associated to PTV_BNCT.'])
                disp('*****')
            end
        end
        
    end
    tissueBin(size(tissueBin,2)+1).name = 'bnct_material';
    tissueBin(size(tissueBin,2)).matIndex = size(tissueBin,2);
    tissueBin(size(tissueBin,2)).linIndVol = lindInd_PTV_BNCT{1};
    [dumIdx1, dumIdx2, dumIdx3] = ind2sub(size(CTdata), lindInd_PTV_BNCT{1});
    tissueBin(size(tissueBin,2)).matIndVol = [dumIdx1, dumIdx2, dumIdx3];
elseif ~exist('bnct_Control')
    disp('*****')
    disp('No PTV for BNCT detected.')
    disp('*****')
    % Add bnct tissue as tissue type
    bnctBinIdx = size(tissueBin,2)+1;
    tissueBin(bnctBinIdx).name = 'bnct_material';
    tissueBin(bnctBinIdx).matIndex = bnctBinIdx;
    tissueBin(bnctBinIdx).linIndVol = [];
    tissueBin(bnctBinIdx).matIndVol = [];
end

%% Appendix: add bone to cst structure
% Generate matRad structure in cst variable
cstCounter = size(cst,1)+1;
cst{cstCounter,1} = cstCounter-1;
cst{cstCounter,2} = 'BONE';
cst{cstCounter,3} = 'OAR';
cst{cstCounter,4} = {tissueBin(boneIndex).linIndVol};
cst{cstCounter,5}.Visible = true;
cst{cstCounter,5}.Priority = 2;

% Make sure linear and sub-indices are the same
for counterTissueBins = 1:size(tissueBin,2)
    tissueBin(counterTissueBins).matIndVol = [];
    if ~isempty(tissueBin(counterTissueBins).linIndVol)
        [tissueBin(counterTissueBins).matIndVol(:,1), tissueBin(counterTissueBins).matIndVol(:,2), tissueBin(counterTissueBins).matIndVol(:,3)] = ind2sub(size(CTdata), tissueBin(counterTissueBins).linIndVol);
    end
end


%% Check if all voxels are still assigned to a medium
dummyVoxelSum = 0;
for checkerC = 1:size(tissueBin,2)
    dummyVoxelSum = dummyVoxelSum + length(tissueBin(checkerC).linIndVol);
end

if dummyVoxelSum > numel(CTdata)
    error('Something went wrong with the segmentation process. Too many voxels after segmentation.')
elseif dummyVoxelSum < numel(CTdata)
    error('Something went wrong with the segmentation process. Voxels were lost in segmentation process.')
else
    disp('Segmentation process performed sucessfully.')
end
