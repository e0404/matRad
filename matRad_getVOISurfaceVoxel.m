function VOI2SurfaceVoxelIDs = matRad_getVOISurfaceVoxel(ct,VOIVoxelIDs)

VOIMask              = zeros(ct.cubeDim);
VOISurfaceMask       = zeros(ct.cubeDim);
VOIMask(VOIVoxelIDs) = 1;

% along x
for i = 1:ct.cubeDim(2)
    
    C = contourc(squeeze(VOIMask(:,i,:)),[1 1]); 
    
    if ~isempty(C)

        % extract coordinates
        logicalIdx = [false,true(1,C(2,1))];
        while length(logicalIdx) < size(C,2)
            logicalIdx = [logicalIdx,false,true(1,C(2,length(logicalIdx) + 1))];
        end
        
        C = C(:,logicalIdx);
        
        % set VOI surface voxels to 1
        for j = 1:size(C,2)
            VOISurfaceMask(C(2,j),i,C(1,j)) = 1;
        end
    end    
end

% along y
for i = 1:ct.cubeDim(1)
    
    C = contourc(squeeze(VOIMask(i,:,:)),[1 1]); 
    
    if ~isempty(C)

        % extract coordinates
        logicalIdx = [false,true(1,C(2,1))];
        while length(logicalIdx) < size(C,2)
            logicalIdx = [logicalIdx,false,true(1,C(2,length(logicalIdx) + 1))];
        end
        
        C = C(:,logicalIdx);
        
        % set VOI surface voxels to 1
        for j = 1:size(C,2)
            VOISurfaceMask(i,C(2,j),C(1,j)) = 1;
        end
    end    
end

% along z
for i = 1:ct.cubeDim(3)
    
    C = contourc(squeeze(VOIMask(:,:,i)),[1 1]); 
    
    if ~isempty(C)

        % extract coordinates
        logicalIdx = [false,true(1,C(2,1))];
        while length(logicalIdx) < size(C,2)
            logicalIdx = [logicalIdx,false,true(1,C(2,length(logicalIdx) + 1))];
        end
        
        C = C(:,logicalIdx);
        
        % set VOI surface voxels to 1
        for j = 1:size(C,2)
            VOISurfaceMask(C(2,j),C(1,j),i) = 1;
        end
    end    
end

VOI2SurfaceVoxelIDs = find(VOISurfaceMask);

end