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
    
    cubeIso = matRad_world2cubeCoords(meanIsoCenter,ct);
    
    % find radius of inner circle from isocenter
    r = 0.8*min([abs([1 ct.cubeDim(1)]-cubeIso(1)) abs([1 ct.cubeDim(2)]-cubeIso(2))]);
    
    % coordinates of circle
    x = r*cosd(0:360)+cubeIso(1);
    y = r*sind(0:360)+cubeIso(2);

    gantryAngleVisColor = 'w';

    hold(axesHandle,'on');
    plot(axesHandle,x,y,'LineWidth',1,'Color',gantryAngleVisColor)

    % add text
    txt = '180°';
    text(axesHandle,1.1*r*sind(0)+cubeIso(1),1.1*r*cosd(0)+cubeIso(2),txt,'Color',gantryAngleVisColor)
    txt = '90°';
    text(axesHandle,1.1*r*sind(90)+cubeIso(1),1.1*r*cosd(90)+cubeIso(2),txt,'Color',gantryAngleVisColor)
    txt = '0°';
    text(axesHandle,1.1*r*sind(180)+cubeIso(1),1.1*r*cosd(180)+cubeIso(2),txt,'Color',gantryAngleVisColor)
    txt = '270°';
    text(axesHandle,1.22*r*sind(270)+cubeIso(1),1.22*r*cosd(270)+cubeIso(2),txt,'Color',gantryAngleVisColor)

    % plot gantry angles
    for i = 1:numel(pln.propStf.gantryAngles)
        plot(axesHandle,[0 r*sind(180-pln.propStf.gantryAngles(i))]+cubeIso(1),[0 r*cosd(180-pln.propStf.gantryAngles(i))]+cubeIso(2),'LineWidth',1,'Color',gantryAngleVisColor)
    end
    
end