function matRad_plotPlan3D(axesHandle,pln,stf)
%MATRAD_PLOTPLAN3D Summary of this function goes here
%   Detailed explanation goes here

vectors = cell(numel(stf),2);

hold(axesHandle,'on');


%draws rays as lines
for fieldIx = 1:numel(stf)
    beamTarget = stf(fieldIx).isoCenter;
    beamSource = stf(fieldIx).sourcePoint + stf(fieldIx).isoCenter;
    beamVector = beamTarget - beamSource;
    
    %{
    for rayIx = 1:numel(stf(fieldIx).ray)
        ray = stf(fieldIx).ray(rayIx);
        rayTarget = ray.targetPoint + stf(fieldIx).isoCenter;
        rayVector = rayTarget - beamSource;
        
        line([beamSource(1) rayTarget(1)],[beamSource(2) rayTarget(2)],[beamSource(3) rayTarget(3)],'Parent',axesHandle,'LineStyle','-','Color',0.5*[1 1 1])
    end
    %}
    
    
    %Alternate approach
    bixelWidth = stf(fieldIx).bixelWidth;
    %Accumulate ray positions in matrix
    rayPos = [];
    rayPos = [stf(fieldIx).ray(:).rayPos_bev];
    rayPos = reshape(rayPos,[3 numel(rayPos)/3])';
    
    %Compute a ray matrix with ones where a ray is
    %Maximum absolute values give extent in one direction
    symmMaxExtent = max([abs(min(rayPos)); abs(max(rayPos))]);
    %we want to have indices and not mm
    symmMaxExtent = symmMaxExtent ./ bixelWidth;
    %now make it symmetric
    symmMaxExtent = 2*symmMaxExtent + [1 0 1];
    %extra padding of one element in each direction to handle contours right
    symmMaxExtent = symmMaxExtent + [2 0 2];
    
    rayMat = zeros(symmMaxExtent(1),symmMaxExtent(3));
    for rayIx = 1:size(rayPos,1)
        centerMat = (size(rayMat)-1) ./ 2 + 1;
        el2D = rayPos(rayIx,[1 3]);
        el2D = el2D ./ bixelWidth;
        el2DIx = centerMat + el2D;
        rayMat(el2DIx(1),el2DIx(2)) = 1;
    end
    %Create the x and y vectors of the ray matrix in world space
    
    %Create contour
    fieldContour2D = contourc(rayMat,1);
    
    %Column in the contour matrix
    cColumn = 1;
    contourIx = 0;
        
    while cColumn <= size(fieldContour2D,2)
        contourIx = contourIx + 1;
        
        %Get contour data
        fieldContour3D(contourIx).level = fieldContour2D(1,cColumn);
        fieldContour3D(contourIx).numVertices = fieldContour2D(2,cColumn);
        fieldContour3D(contourIx).coord_bev = [ fieldContour2D(2,cColumn+1:cColumn+fieldContour3D(contourIx).numVertices); ... x
                                                zeros(1,fieldContour3D(contourIx).numVertices);
                                                fieldContour2D(1,cColumn+1:cColumn+fieldContour3D(contourIx).numVertices)]; ... y];
                                                
        fieldContour3D(contourIx).coord_bev = bsxfun(@minus,fieldContour3D(contourIx).coord_bev,symmMaxExtent'./2) * bixelWidth;
        
        % Transform to world space
        % compute coordinates in lps coordinate system, i.e. rotate beam
        % geometry around fixed patient; 
    
        % Rotation around Z axis (gantry)
        rotMx_XY = [    cosd(pln.gantryAngles(fieldIx))   -sind(pln.gantryAngles(fieldIx))  0;
                        sind(pln.gantryAngles(fieldIx))    cosd(pln.gantryAngles(fieldIx))  0;
                        0                            0                          1];
    
        % Rotation around Y axis (couch)
        rotMx_XZ = [    cosd(pln.couchAngles(fieldIx))    0    sind(pln.couchAngles(fieldIx));
                        0                           1    0;
                       -sind(pln.couchAngles(fieldIx))    0    cosd(pln.couchAngles(fieldIx))];
        
        rotMat = rotMx_XY *rotMx_XZ;
                   
        for v = 1:fieldContour3D(contourIx).numVertices            
            fieldContour3D(contourIx).coord(1:3,v) = rotMat * fieldContour3D(contourIx).coord_bev(1:3,v) + beamTarget';
            line('XData',[beamSource(1) fieldContour3D(contourIx).coord(1,v)],'YData',[beamSource(2) fieldContour3D(contourIx).coord(2,v)],'ZData',[beamSource(3) fieldContour3D(contourIx).coord(3,v)],'Parent',axesHandle,'LineWidth',1,'Color',[255 20 147]/255);
        end    
        
        %Draw the contour in the isocenter
        line('XData',fieldContour3D(contourIx).coord(1,:),'YData',fieldContour3D(contourIx).coord(2,:),'ZData',fieldContour3D(contourIx).coord(3,:),'Parent',axesHandle,'LineWidth',2,'Color',[255 20 147]/255);
        
        %Draw the lines from the source point
        
        
        cColumn = cColumn + fieldContour3D(contourIx).numVertices + 1; 
    end
    
end


%{
for fieldIx = 1:numel(stf)
    rayPos_bev = cat(1,stf(fieldIx).ray.rayPos_bev);
    [~,minPosIx] = min(rayPos_bev);
    [~,maxPosIx] = max(rayPos_bev);
    
    minPosX = stf(fieldIx).ray(minPosIx(1)).rayPos + stf(fieldIx).isoCenter;
    minPosZ = stf(fieldIx).ray(minPosIx(3)).rayPos + stf(fieldIx).isoCenter;
    maxPosX = stf(fieldIx).ray(maxPosIx(1)).rayPos + stf(fieldIx).isoCenter;
    maxPosZ = stf(fieldIx).ray(maxPosIx(3)).rayPos + stf(fieldIx).isoCenter;
    
    planeX = [minPosX(1) minPosZ(1) maxPosZ(1) maxPosX(1)];
    planeY = [minPosX(2) minPosZ(2) maxPosZ(2) maxPosX(2)];
    planeZ = [minPosX(3) minPosZ(3) maxPosZ(3) maxPosX(3)];
    
    fill3(axesHandle,planeX,planeY,planeZ,'red');
end
%}
end



