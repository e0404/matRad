function s = matRad_shadowPlot(x, yLow, yUp, color, legendName, alphaTrans)
% shadowPlot to plot confidence bands mainly
% 
% call
%   shadowPlot(x, yLow, yUp, color, legendName, alphaTrans)
%
% input
%   x:          x axis values
%   yLow:       lower bound (start of shadowing)
%   yUp:        upper bound (end of shadowing)
%   color:      color as [R G B]
%   legendName: legendname to be shown in the plot
%   alphaTrans: transparency
%
% output
%   axesHandle
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set alpha value to a default if not parsed
if ~exist('alphaTrans', 'var') || isempty(alphaTrans)
    alphaTrans = .3;
end

% x, yLow, yUp has to be col vector
if size(x, 1) == 1
    x = x';
end
if size(yLow, 1) == 1
    yLow = yLow';
end
if size(yUp, 1) == 1
    yUp = yUp';
end

% check for nan
if any(isnan(yLow)) || any(isnan(yUp))
    s = gca;
    return
end

% no shadow under the lower limit line
sThick.low = 0 * yLow;
sThick.up = yUp - yLow;    

s = fill([x;flipud(x)],[yLow - sThick.low;flipud(yLow + sThick.up)], color,'linestyle','none', 'DisplayName', legendName);
alpha(s, alphaTrans)
hold on;

end % eof
