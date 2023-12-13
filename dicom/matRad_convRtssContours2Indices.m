function indices = matRad_convRtssContours2Indices(structure,ct)
% matRad function to convert a polygon segmentation from an rt structure
% set into a binary segmentation as required within matRad's cst struct
% 
% call
%   indices = matRad_convRtssContours2Indices(contPoints,ct)
%
% input
%   structure:      information about a single structure
%   ct:             matRad ct struct where the binary segmentations will
%                   be aligned to
%
% output
%   indicies:       indices of voxels of the ct cube that are inside the
%                   contour
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
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
voiCube = zeros(ct.cubeDim);

% loop over all closed contour items
for i = 1:size(structure.item,2)

    if ~isempty(structure.item(i).points)

        dicomCtSlicePos = unique(structure.item(i).points(:,3));
        
        if numel(dicomCtSlicePos) > 1 || isempty(dicomCtSlicePos)
            matRad_cfg.dispError('Contour defined over multiple planes!');
        end
    
        round2 = @(a,b) round(a*10^b)/10^b;
        dicomCtSliceThickness = ct.dicomInfo.SliceThickness(round2(ct.dicomInfo.SlicePositions,1)==round2(dicomCtSlicePos,1));
        
        %Sanity check
        msg = checkSliceThickness(dicomCtSliceThickness);
        if ~isempty(msg)
            matRad_cfg.dispError('Slice Thickness of slice at %f could not be identified: %s',dicomCtSlicePos,msg);
        end
        
        slicesInMatradCt = find(dicomCtSlicePos+dicomCtSliceThickness/2 > ct.z & dicomCtSlicePos-dicomCtSliceThickness/2 <= ct.z);
        
        coords1 = interp1(ct.x,1:ct.cubeDim(2),structure.item(i).points(:,1),'linear','extrap');
        coords2 = interp1(ct.y,1:ct.cubeDim(1),structure.item(i).points(:,2),'linear','extrap');
        
        binIn = poly2mask(coords1,coords2,ct.cubeDim(1),ct.cubeDim(2));
        
        % loop over all slices in matRad ct
        for j = 1:numel(slicesInMatradCt)
            voiCube(:,:,slicesInMatradCt(j)) = voiCube(:,:,slicesInMatradCt(j)) | binIn;
        end
        
    end
    
end
% The x- & y-direction in lps-coordinates are specified in:
% ImageOrientationPatient

xDir = ct.dicomInfo.ImageOrientationPatient(1:3); % lps: [1;0;0]
yDir = ct.dicomInfo.ImageOrientationPatient(4:6); % lps: [0;1;0]
nonStandardDirection = false;

if xDir(1) == 1 && xDir(2) == 0 && xDir(3) == 0
%     matRad_cfg.dispInfo('x-direction OK\n')
elseif xDir(1) == -1 && xDir(2) == 0 && xDir(3) == 0
    voiCube = flip(voiCube,1);
else
    nonStandardDirection = true;
end
    
if yDir(1) == 0 && yDir(2) == 1 && yDir(3) == 0
%     matRad_cfg.dispInfo('y-direction OK\n')
elseif yDir(1) == 0 && yDir(2) == -1 && yDir(3) == 0
    voiCube = flip(voiCube,2);
else
    nonStandardDirection = true;
end

if nonStandardDirection
    matRad_cfg.dispWarning(['Non-standard patient orientation.\n'...
        'CT might not fit to contoured structures\n'])
end

indices = find(voiCube(:));

end

function msg = checkSliceThickness(dicomCtSliceThickness)
    if isempty(dicomCtSliceThickness)
        msg = 'Slice could not be identified (empty)';
    elseif ~isscalar(dicomCtSliceThickness)
        msg = 'Slice thickness not unique';
    elseif ~isnumeric(dicomCtSliceThickness)
        msg = 'unexpected value';
    else
        msg = '';
    end
end
