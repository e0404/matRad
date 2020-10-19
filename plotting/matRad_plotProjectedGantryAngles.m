function matRad_plotProjectedGantryAngles(axesHandle,pln,ct,plane)
% matRad function that plots all gantry angles 
% projected to the coplanar plane if current view is axial view
%
% call
%   matRad_plotProjectedGantryAngles(axesHandle,pln,ct,plane)
%
% input
%   axesHandle: handle to the axis where the plot shouldd appear
%   pln:        matRad pln struct
%   ct:         matRad ct struct
%   plane:      current view plane
%
% output
%   -
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% show only in axial plane
if plane == 3
   
    meanIsoCenter = mean(pln.propStf.isoCenter,1);
    
    xOffset = meanIsoCenter(1)/ct.resolution.x;
    yOffset = meanIsoCenter(2)/ct.resolution.y;
    
    % find radius of inner circle from isocenter
    r = 0.8*min([abs([1 ct.cubeDim(1)]-xOffset) abs([1 ct.cubeDim(2)]-yOffset)]);
    
    % coordinates of circle
    x = r*cosd(0:360)+xOffset;
    y = r*sind(0:360)+yOffset;

    gantryAngleVisColor = 'w';

    plot(axesHandle,x,y,'LineWidth',1,'Color',gantryAngleVisColor)

    % add text
    txt = '180°';
    text(axesHandle,1.1*r*sind(0)+xOffset,1.1*r*cosd(0)+yOffset,txt,'Color',gantryAngleVisColor)
    txt = '90°';
    text(axesHandle,1.1*r*sind(90)+xOffset,1.1*r*cosd(90)+yOffset,txt,'Color',gantryAngleVisColor)
    txt = '0°';
    text(axesHandle,1.1*r*sind(180)+xOffset,1.1*r*cosd(180)+yOffset,txt,'Color',gantryAngleVisColor)
    txt = '270°';
    text(axesHandle,1.22*r*sind(270)+xOffset,1.22*r*cosd(270)+yOffset,txt,'Color',gantryAngleVisColor)

    % plot gantry angles
    for i = 1:numel(pln.propStf.gantryAngles)
        plot(axesHandle,[0 r*sind(180-pln.propStf.gantryAngles(i))]+xOffset,[0 r*cosd(180-pln.propStf.gantryAngles(i))]+yOffset,'LineWidth',1,'Color',gantryAngleVisColor)
    end
    
end