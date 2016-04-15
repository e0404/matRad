function VOISurfaceVoxel = matRad_getVOISurfaceVoxel(ct,VOIVoxelIDs)

    % get cube coordinates
    [xCoords_vox, yCoords_vox, zCoords_vox] = meshgrid(1:ct.cubeDim(1),1:ct.cubeDim(2),1:ct.cubeDim(3));

    xCoords = xCoords_vox(:);
    yCoords = yCoords_vox(:);
    zCoords = zCoords_vox(:);

    % get VOI coordinates
    xCoordsVOI = xCoords(VOIVoxelIDs);
    yCoordsVOI = yCoords(VOIVoxelIDs);
    zCoordsVOI = zCoords(VOIVoxelIDs);

    xCoordsVOIUnique = unique(xCoordsVOI); 
    yCoordsVOIUnique = unique(yCoordsVOI); 
    zCoordsVOIUnique = unique(zCoordsVOI); 

    % x
    VOISurfaceX = [];
    for i = 1:length(xCoordsVOIUnique)

        idx = find(xCoordsVOI == xCoordsVOIUnique(i));

        if i ~= 1 && i ~= length(xCoordsVOIUnique)
            idx = xCoordsVOI == xCoordsVOIUnique(i) & yCoordsVOI == min(yCoordsVOI(idx)) |...
                  xCoordsVOI == xCoordsVOIUnique(i) & yCoordsVOI == max(yCoordsVOI(idx)) |...
                  xCoordsVOI == xCoordsVOIUnique(i) & zCoordsVOI == min(zCoordsVOI(idx)) |...
                  xCoordsVOI == xCoordsVOIUnique(i) & zCoordsVOI == max(zCoordsVOI(idx));

            idx = find(idx);  
        end

        VOISurfaceX = [VOISurfaceX;VOIVoxelIDs(idx)];   
    end

    % y
    VOISurfaceY = [];
    for i = 1:length(yCoordsVOIUnique)

        idx = find(yCoordsVOI == yCoordsVOIUnique(i));

        if i ~= 1 && i ~= length(yCoordsVOIUnique)
            idx = yCoordsVOI == yCoordsVOIUnique(i) & xCoordsVOI == min(xCoordsVOI(idx)) |...
                  yCoordsVOI == yCoordsVOIUnique(i) & xCoordsVOI == max(xCoordsVOI(idx)) |...
                  yCoordsVOI == yCoordsVOIUnique(i) & zCoordsVOI == min(zCoordsVOI(idx)) |...
                  yCoordsVOI == yCoordsVOIUnique(i) & zCoordsVOI == max(zCoordsVOI(idx));

            idx = find(idx);  
        end

        VOISurfaceY = [VOISurfaceY;VOIVoxelIDs(idx)];   
    end

    % z
    VOISurfaceZ = [];
    for i = 1:length(zCoordsVOIUnique)

        idx = find(zCoordsVOI == zCoordsVOIUnique(i));

        if i ~= 1 && i ~= length(zCoordsVOIUnique)
            idx = zCoordsVOI == zCoordsVOIUnique(i) & yCoordsVOI == min(yCoordsVOI(idx)) |...
                  zCoordsVOI == zCoordsVOIUnique(i) & yCoordsVOI == max(yCoordsVOI(idx)) |...
                  zCoordsVOI == zCoordsVOIUnique(i) & xCoordsVOI == min(xCoordsVOI(idx)) |...
                  zCoordsVOI == zCoordsVOIUnique(i) & xCoordsVOI == max(xCoordsVOI(idx));

            idx = find(idx);  
        end

        VOISurfaceZ = [VOISurfaceZ;VOIVoxelIDs(idx)];   
    end

    VOISurfaceVoxel = union(VOISurfaceX,VOISurfaceY);
    VOISurfaceVoxel = union(VOISurfaceVoxel,VOISurfaceZ); 

end