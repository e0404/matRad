function indices = matRad_convRtssContours2Indices(contPoints,ct)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to convert a polygon segmentation from an rt structure
% set into a binary segmentation as required within matRad's cst struct
% 
% call
%   indices = matRad_convRtssContours2Indices(contPoints,ct)
%
% input
%   contPoints:     set of all contour points belonging to a single
%                   structure
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

slicesInOrigDicomWithContours = unique(contPoints(:,3));

[X,Y] = meshgrid(ct.x,ct.y);

% loop over all slices where contour points have been defined in the original dicom ct
for i = 1:numel(slicesInOrigDicomWithContours)

    dicomCtSlicePos       = slicesInOrigDicomWithContours(i);
    
    dicomCtSliceThickness = ct.dicomInfo.SliceThickness(round(ct.dicomInfo.SlicePositions,2)==dicomCtSlicePos);
    
    contPointsInCurrDicomSlice = contPoints(contPoints(:,3) == dicomCtSlicePos,:);
    
    binIn = inpolygon(X,Y,contPointsInCurrDicomSlice(:,1),contPointsInCurrDicomSlice(:,2));
    
    slicesInMatradCt = find(dicomCtSlicePos+dicomCtSliceThickness/2 > ct.z & dicomCtSlicePos-dicomCtSliceThickness/2 <= ct.z);
        
    % loop over all slices in matRad ct
    for j = 1:numel(slicesInMatradCt)
        voiCube(:,:,slicesInMatradCt(j)) = binIn;
    end

end

indices = find(voiCube>0);
