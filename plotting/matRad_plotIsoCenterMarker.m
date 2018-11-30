function markerHandle = matRad_plotIsoCenterMarker(axesHandle,pln,ct,plane,slice,style)
% matRad function that plots an isocenter marker
%
% call
%   markerHandle = matRad_plotIsoCenterMarker(axesHandle,pln,ct,plane,slice)
%
% input
%   axesHandle          handle to axes the marker should be displayed in
%   pln                 matRad plan structure
%   ct                  matRad ct structure
%   plane               plane view (coronal=1,sagittal=2,axial=3)
%   slice               slice in the selected plane of the 3D cube
%   style               optional argument can be 'marker' or 'lines'. 
%                       this might be useful for some issues 
%                       with export to tikz for example
%                       
%
% output
%   markerHandle:       handle to the isocenter marker
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

uniqueIsoCenters = unique(pln.propStf.isoCenter,'rows');

for i = 1:size(uniqueIsoCenters,1)

    vIsoCenter           = round(uniqueIsoCenters(i,:)./[ct.resolution.x ct.resolution.y ct.resolution.z]);
    if  plane == 3% Axial plane
        vIsoCenterPlot  = [vIsoCenter(1) vIsoCenter(2)];
        if vIsoCenter(3) == slice
            isoCenterDirection = 0;
        else 
            isoCenterDirection = sign(double(vIsoCenter(3) > slice) - 0.5);
        end
    elseif plane == 2
        vIsoCenterPlot  = [vIsoCenter(3) vIsoCenter(2)];
        if vIsoCenter(2) == slice
            isoCenterDirection = 0;
        else
            isoCenterDirection = sign(double(vIsoCenter(1) > slice) - 0.5);
        end

    elseif plane == 1    
        vIsoCenterPlot  = [vIsoCenter(3) vIsoCenter(1)];
        if vIsoCenter(1) == slice
            isoCenterDirection = 0;
        else
            isoCenterDirection = sign(double(vIsoCenter(2) > slice) - 0.5);
        end
    end



    if nargin < 6
        style = 'lines';
    end

    if strcmpi(style,'lines') == 0 && strcmpi(style,'marker') == 0
        warning('Style option not recognized. Using lines as default');
        style = 'lines';
    end

    markerSize = 36/ct.resolution.x;
    markerColor = [0.27 0.27 0.27];

    if strcmpi(style,'marker')
        if isoCenterDirection > 0
            markerStyle = '^';
        elseif isoCenterDirection < 0
            markerStyle = 'v';
        else
            markerStyle = 'x';
            markerColor = [0 0 0];
        end    
        markerHandle = plot(axesHandle,vIsoCenterPlot(1),vIsoCenterPlot(2),markerStyle,'MarkerSize',markerSize,'LineWidth',4,'Color',markerColor);
    else
        if isoCenterDirection > 0
            linesX = [vIsoCenterPlot(1)-markerSize/4 vIsoCenterPlot(1) vIsoCenterPlot(1)+markerSize/4];
            linesY = [vIsoCenterPlot(2)+markerSize/4 vIsoCenterPlot(2) vIsoCenterPlot(2)+markerSize/4];
        elseif isoCenterDirection < 0        
            linesX = [vIsoCenterPlot(1)-markerSize/4 vIsoCenterPlot(1) vIsoCenterPlot(1)+markerSize/4];
            linesY = [vIsoCenterPlot(2)-markerSize/4 vIsoCenterPlot(2) vIsoCenterPlot(2)-markerSize/4];
        else
            linesX = [vIsoCenterPlot(1)-markerSize/4 vIsoCenterPlot(1)-markerSize/4; vIsoCenterPlot(1)+markerSize/4 vIsoCenterPlot(1)+markerSize/4];
            linesY = [vIsoCenterPlot(2)-markerSize/4 vIsoCenterPlot(2)+markerSize/4; vIsoCenterPlot(2)+markerSize/4 vIsoCenterPlot(2)-markerSize/4];
            markerColor = [0 0 0];
        end

        markerHandle = line(linesX,linesY,'LineWidth',4,'Color',markerColor,'Parent',axesHandle);
    end

end
