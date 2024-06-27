function colorMap = diffMap(cMapSize)
% matRad difference colormap 
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

colorMapData = ...
    [0         0    1.0000
    0.0323    0.0323    1.0000
    0.0645    0.0645    1.0000
    0.0968    0.0968    1.0000
    0.1290    0.1290    1.0000
    0.1613    0.1613    1.0000
    0.1935    0.1935    1.0000
    0.2258    0.2258    1.0000
    0.2581    0.2581    1.0000
    0.2903    0.2903    1.0000
    0.3226    0.3226    1.0000
    0.3548    0.3548    1.0000
    0.3871    0.3871    1.0000
    0.4194    0.4194    1.0000
    0.4516    0.4516    1.0000
    0.4839    0.4839    1.0000
    0.5161    0.5161    1.0000
    0.5484    0.5484    1.0000
    0.5806    0.5806    1.0000
    0.6129    0.6129    1.0000
    0.6452    0.6452    1.0000
    0.6774    0.6774    1.0000
    0.7097    0.7097    1.0000
    0.7419    0.7419    1.0000
    0.7742    0.7742    1.0000
    0.8065    0.8065    1.0000
    0.8387    0.8387    1.0000
    0.8710    0.8710    1.0000
    0.9032    0.9032    1.0000
    0.9355    0.9355    1.0000
    0.9677    0.9677    1.0000
    1.0000    1.0000    1.0000
    1.0000    1.0000    1.0000
    1.0000    0.9677    0.9677
    1.0000    0.9355    0.9355
    1.0000    0.9032    0.9032
    1.0000    0.8710    0.8710
    1.0000    0.8387    0.8387
    1.0000    0.8065    0.8065
    1.0000    0.7742    0.7742
    1.0000    0.7419    0.7419
    1.0000    0.7097    0.7097
    1.0000    0.6774    0.6774
    1.0000    0.6452    0.6452
    1.0000    0.6129    0.6129
    1.0000    0.5806    0.5806
    1.0000    0.5484    0.5484
    1.0000    0.5161    0.5161
    1.0000    0.4839    0.4839
    1.0000    0.4516    0.4516
    1.0000    0.4194    0.4194
    1.0000    0.3871    0.3871
    1.0000    0.3548    0.3548
    1.0000    0.3226    0.3226
    1.0000    0.2903    0.2903
    1.0000    0.2581    0.2581
    1.0000    0.2258    0.2258
    1.0000    0.1935    0.1935
    1.0000    0.1613    0.1613
    1.0000    0.1290    0.1290
    1.0000    0.0968    0.0968
    1.0000    0.0645    0.0645
    1.0000    0.0323    0.0323
    1.0000         0         0];

if nargin < 1
    colorMap = colorMapData;
elseif size(colorMapData,1) == cMapSize
    colorMap = colorMapData;
else
    %We have to interpolate the colormap
    newX = linspace(1,64,cMapSize);
    oldX = 1:64;
    colorMap = interp1(oldX,colorMapData,newX);
    %{
    %resample via HSV.. more color-true than above, but doesn't work with
    %every colormap
    hsv                        = rgb2hsv(cm);
    hsv(144:end,1)             = hsv(144:end,1)+1;
    ColorMap                   = interp1(linspace(0,1,size(cm,1)),hsv,linspace(0,1,cMapSize));
    ColorMap(cm(:,1)>1,1) = ColorMap(cm(:,1)>1,1)-1;
    ColorMap                   = hsv2rgb(cm);
    %}
end

end

