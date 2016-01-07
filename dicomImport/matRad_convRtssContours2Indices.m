function indices = matRad_convRtssContours2Indices(structure,ct)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

voiCube = zeros(size(ct.cube));

[X,Y] = meshgrid(ct.x,ct.y);

% loop over all closed contour items
for i = 1:size(structure.item,2)

    if ~isempty(structure.item(i).points)

        dicomCtSlicePos = unique(structure.item(i).points(:,3));
        
        if numel(dicomCtSlicePos) > 1
            error('Contour defined over multiple planes\n');
        end
    
        round2 = @(a,b) round(a*10^b)/10^b;
        dicomCtSliceThickness = ct.dicomInfo.SliceThickness(round2(ct.dicomInfo.SlicePositions,2)==round2(dicomCtSlicePos,2));

        binIn = inpolygon(X,Y,structure.item(i).points(:,1),structure.item(i).points(:,2));

        slicesInMatradCt = find(dicomCtSlicePos+dicomCtSliceThickness/2 > ct.z & dicomCtSlicePos-dicomCtSliceThickness/2 <= ct.z);

        % loop over all slices in matRad ct
        for j = 1:numel(slicesInMatradCt)
            voiCube(:,:,slicesInMatradCt(j)) = voiCube(:,:,slicesInMatradCt(j)) | binIn;
        end
        
    end
    
end

indices = find(voiCube>0);
