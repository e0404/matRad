function stf = matRad_shiftStf(stf,ct,offset)

% check input data
if size(offset,1) ~= 1 && size(offset,2) ~= 3
    error('offset must be a 1x3 row vector\n');
end

DensityThresholdSSD = 0.05;

for i = 1:size(stf,2)
    
    % rotation around Z axis (gantry)
    inv_rotMx_XY_T = [ cosd(-stf(i).gantryAngle) sind(-stf(i).gantryAngle) 0;
                      -sind(-stf(i).gantryAngle) cosd(-stf(i).gantryAngle) 0;
                                                  0                            0 1];

    % rotation around Y axis (couch)
    inv_rotMx_XZ_T = [cosd(-stf(i).couchAngle) 0 -sind(-stf(i).couchAngle);
                                                0 1                            0;
                      sind(-stf(i).couchAngle) 0  cosd(-stf(i).couchAngle)];

    % Rotate offset into bev
    offset_bev = offset*inv_rotMx_XZ_T*inv_rotMx_XY_T;
    
    % add offset
    stf(i).isoCenter       = stf(i).isoCenter + offset;
    stf(i).sourcePoint     = stf(i).sourcePoint + offset;
    stf(i).sourcePoint_bev = stf(i).sourcePoint_bev + offset_bev;
    
    for j = 1:stf(i).numOfRays
        
        stf(i).ray(j).rayPos_bev      = stf(i).ray(j).rayPos_bev + offset_bev;
        stf(i).ray(j).targetPoint_bev = stf(i).ray(j).targetPoint_bev + offset_bev;
        stf(i).ray(j).rayPos          = stf(i).ray(j).rayPos + offset;
        stf(i).ray(j).targetPoint     = stf(i).ray(j).targetPoint + offset;
        
        % ray tracing necessary to determine depth of the target
        [alpha,~,rho,~,~] = matRad_siddonRayTracer(stf(i).isoCenter, ...
                             ct.resolution, ...
                             stf(i).sourcePoint, ...
                             stf(i).ray(j).targetPoint, ...
                             {ct.cube});

        ixSSD = find(rho{1} > DensityThresholdSSD,1,'first');

        if isempty(ixSSD)== 1
            warning('Surface for SSD calculation starts directly in first voxel of CT\n');
        end
        
        % calculate SSD
        stf(i).ray(j).SSD = 2 * stf(i).SAD * alpha(ixSSD);
        
    end

end
