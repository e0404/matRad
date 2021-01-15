function matRad_plotPlan3D(axesHandle,pln,stf)
% matRad function to visualize a plan in 3D 
% Stf is optional for plotting more detailed field contours in 
% visualization of the impinging beams.
% 
% call
%  rotMat = matRad_plotPlan3D(axesHandle,pln)
%  rotMat = matRad_plotPlan3D(axesHandle,pln,stf)
%
% input
%   axesHandle: handle to the axes the plan should be visualized in.
%   pln:        matRad plan meta information struct
%   stf:        optional steering information struct. if stf is passed and 
%               not empty, the function will use the ray position  
%               information to plot more detailed field contours than with 
%               pln only
%
% output
%   -
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

if ishold(axesHandle)
    wasHold = true;
else
    wasHold = false;
    hold(axesHandle,'on');
end

%nice pink ;)
beamColor = [255 20 147]/255;

%We perform a rudimentary visualization of the beam angles if there is no stf
if nargin < 3 || isempty(stf)
    visFieldExtent = 50;
    visFieldPoints_bev = [-visFieldExtent 0 -visFieldExtent; ...
        visFieldExtent  0 -visFieldExtent; ...
        visFieldExtent  0  visFieldExtent; ...
        -visFieldExtent 0  visFieldExtent;
        -visFieldExtent 0 -visFieldExtent]; %repeat the first value for closed contour
    
    visFieldPoints_bev = visFieldPoints_bev ./ 2;
    visFieldPoints_bev = visFieldPoints_bev';
    
    fileName = [pln.radiationMode '_' pln.machine];
    %Get a SAD
    try
        load([fileparts(mfilename('fullpath')) filesep 'basedata' filesep fileName]);
        SAD = machine.meta.SAD;
    catch
        if strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
            SAD = 10000;
        else
            SAD = 1500;
        end
    end
    
    % default beam vector
    beamVector = [0 SAD 0];
   
    
    for beamIx = 1:pln.propStf.numOfBeams
        rotMat = matRad_getRotationMatrix(pln.propStf.gantryAngles(beamIx),pln.propStf.couchAngles(beamIx));
        beamIsoCenter = pln.propStf.isoCenter(beamIx,:);
        currBeamVector = rotMat*beamVector';        
        currBeamSource = beamIsoCenter - currBeamVector';
        currBeamOuterTarget = beamIsoCenter + currBeamVector';
        
        %Central ray
        line('XData',[currBeamSource(1) beamIsoCenter(1)],'YData',[currBeamSource(2) beamIsoCenter(2)],'ZData',[currBeamSource(3) beamIsoCenter(3)],'Parent',axesHandle,'LineWidth',2,'Color',beamColor);
        line('XData',[currBeamOuterTarget(1) beamIsoCenter(1)],'YData',[currBeamOuterTarget(2) beamIsoCenter(2)],'ZData',[currBeamOuterTarget(3) beamIsoCenter(3)],'Parent',axesHandle,'LineWidth',2,'Color',beamColor,'LineStyle',':');
        
        %Field rays
        for v = 1:size(visFieldPoints_bev,2)
            visFieldPoints_world(1:3,v) = rotMat*visFieldPoints_bev(1:3,v) + beamIsoCenter';
            
            %skipt the drawing of the first vertex since it is there twice
            if v==1 
                continue;
            end
            
            
            
            %Draw the lines from the source point to the contour point           
            line('XData',[currBeamSource(1) visFieldPoints_world(1,v)],'YData',[currBeamSource(2) visFieldPoints_world(2,v)],'ZData',[currBeamSource(3) visFieldPoints_world(3,v)],'Parent',axesHandle,'LineWidth',1,'Color',beamColor);
            %extend further
            currPointVector = visFieldPoints_world(1:3,v) - currBeamSource';
            currPointOuterTarget = currPointVector + visFieldPoints_world(1:3,v);
            line('XData',[currPointOuterTarget(1) visFieldPoints_world(1,v)],'YData',[currPointOuterTarget(2) visFieldPoints_world(2,v)],'ZData',[currPointOuterTarget(3) visFieldPoints_world(3,v)],'Parent',axesHandle,'LineWidth',1,'Color',beamColor,'LineStyle',':');
            
            %contour
            line('XData',[visFieldPoints_world(1,v) visFieldPoints_world(1,v-1)],...
                'YData',[visFieldPoints_world(2,v) visFieldPoints_world(2,v-1)],...
                'ZData',[visFieldPoints_world(3,v) visFieldPoints_world(3,v-1)],...
                'Parent',axesHandle,'LineWidth',2,'Color',beamColor);
        end
    end
else %We use the steering information to visualize the field contour
    for fieldIx = 1:numel(stf)
        beamTarget = stf(fieldIx).isoCenter;
        beamSource = stf(fieldIx).sourcePoint + stf(fieldIx).isoCenter;
        
        rotMat = matRad_getRotationMatrix(pln.propStf.gantryAngles(fieldIx),pln.propStf.couchAngles(fieldIx));
        
        bixelWidth = stf(fieldIx).bixelWidth;
        %Accumulate ray positions in matrix
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
        
        %Fill the ray matrix with ones at positions of rays
        rayMat = zeros(symmMaxExtent(1),symmMaxExtent(3));
        for rayIx = 1:size(rayPos,1)
            centerMat = (size(rayMat)-1) ./ 2 + 1;
            el2D = rayPos(rayIx,[1 3]);
            el2D = el2D ./ bixelWidth;
            el2DIx = centerMat + el2D;
            rayMat(el2DIx(1),el2DIx(2)) = 1;
        end
        %Create contour of the field
        fieldContour2D = contourc(rayMat,1);
        
        %Column in the contour matrix
        cColumn = 1;
        contourIx = 0;
        
        %Orientate the contour in 3D world space
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
            
            for v = 1:fieldContour3D(contourIx).numVertices
                fieldContour3D(contourIx).coord(1:3,v) = rotMat * fieldContour3D(contourIx).coord_bev(1:3,v) + beamTarget';
                %Draw the lines from the source point to the contour point
                line('XData',[beamSource(1) fieldContour3D(contourIx).coord(1,v)],'YData',[beamSource(2) fieldContour3D(contourIx).coord(2,v)],'ZData',[beamSource(3) fieldContour3D(contourIx).coord(3,v)],'Parent',axesHandle,'LineWidth',1,'Color',beamColor);
            end
            
            %Draw the contour in the isocenter (shows the field)
            line('XData',fieldContour3D(contourIx).coord(1,:),'YData',fieldContour3D(contourIx).coord(2,:),'ZData',fieldContour3D(contourIx).coord(3,:),'Parent',axesHandle,'LineWidth',2,'Color',beamColor);
            
            cColumn = cColumn + fieldContour3D(contourIx).numVertices + 1;
        end
        
    end
end

if ~wasHold
    hold(axesHandle,'off');
end

end
