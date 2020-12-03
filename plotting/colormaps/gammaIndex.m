function colorMap = gammaIndex(cMapSize)
% matRad gamma index colormap  
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
colorMapData = [  0.2081    0.1663    0.5292
                    0.2336    0.1932    0.5444
                    0.2592    0.2201    0.5596
                    0.2847    0.2470    0.5748
                    0.3103    0.2739    0.5899
                    0.3358    0.3008    0.6051
                    0.3614    0.3277    0.6203
                    0.3869    0.3546    0.6355
                    0.4125    0.3814    0.6507
                    0.4380    0.4083    0.6659
                    0.4636    0.4352    0.6811
                    0.4891    0.4621    0.6963
                    0.5146    0.4890    0.7114
                    0.5402    0.5159    0.7266
                    0.5657    0.5428    0.7418
                    0.5913    0.5697    0.7570
                    0.6168    0.5966    0.7722
                    0.6424    0.6235    0.7874
                    0.6679    0.6504    0.8026
                    0.6935    0.6773    0.8178
                    0.7190    0.7042    0.8329
                    0.7445    0.7311    0.8481
                    0.7701    0.7580    0.8633
                    0.7956    0.7849    0.8785
                    0.8212    0.8117    0.8937
                    0.8467    0.8386    0.9089
                    0.8723    0.8655    0.9241
                    0.8978    0.8924    0.9393
                    0.9234    0.9193    0.9544
                    0.9489    0.9462    0.9696
                    0.9745    0.9731    0.9848
                    1.0000    1.0000    1.0000
                    1.0000    1.0000         0
                    1.0000    0.9677         0
                    1.0000    0.9355         0
                    1.0000    0.9032         0
                    1.0000    0.8710         0
                    1.0000    0.8387         0
                    1.0000    0.8065         0
                    1.0000    0.7742         0
                    1.0000    0.7419         0
                    1.0000    0.7097         0
                    1.0000    0.6774         0
                    1.0000    0.6452         0
                    1.0000    0.6129         0
                    1.0000    0.5806         0
                    1.0000    0.5484         0
                    1.0000    0.5161         0
                    1.0000    0.4839         0
                    1.0000    0.4516         0
                    1.0000    0.4194         0
                    1.0000    0.3871         0
                    1.0000    0.3548         0
                    1.0000    0.3226         0
                    1.0000    0.2903         0
                    1.0000    0.2581         0
                    1.0000    0.2258         0
                    1.0000    0.1935         0
                    1.0000    0.1613         0
                    1.0000    0.1290         0
                    1.0000    0.0968         0
                    1.0000    0.0645         0
                    1.0000    0.0323         0
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

