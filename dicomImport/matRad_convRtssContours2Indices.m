function indices = matRad_convRtssContours2Indices(contPoints,ct)

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
