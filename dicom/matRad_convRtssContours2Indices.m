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

voiCube = zeros(ct.cubeDim);

% loop over all closed contour items
for i = 1:size(structure.item,2)

    if ~isempty(structure.item(i).points)

        dicomCtSlicePos = unique(structure.item(i).points(:,3));
        
        if numel(dicomCtSlicePos) > 1 || isempty(dicomCtSlicePos)
            error('Contour defined over multiple planes!');
        end
    
        round2 = @(a,b) round(a*10^b)/10^b;
        dicomCtSliceThickness = ct.dicomInfo.SliceThickness(round2(ct.dicomInfo.SlicePositions,1)==round2(dicomCtSlicePos,1));
        
        %Sanity check
        msg = checkSliceThickness(dicomCtSliceThickness);
        if ~isempty(msg)
            error('Slice Thickness of slice at %f could not be identified: %s',dicomCtSlicePos,msg);
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
