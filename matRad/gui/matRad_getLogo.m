function [im,alpha] = matRad_getLogo(scale)
% matRad function to obtain matRad logo image data adhering to theme
% 
% call
%   matRad_getLogoDKFZ()
%   matRad_getLogoDKFZ(scale)
%
% input:
%   scale:  (optional, default 1): either scalar resizing factor or new 
%           image size that should be used in pixels
%
% output:
%           im:     RGB image
%           alpha:  alpha channel
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2024 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 1
    scale = 1;
end

matRad_cfg = MatRad_Config.instance();
filepath = matRad_cfg.matRadSrcRoot;

%Read Logo
[im,~,alpha] = imread(fullfile(filepath,'gfx','matrad_logo.png'));

%Recolor Logo
if matRad_cfg.isOctave
    rgbUint = uint8(matRad_cfg.gui.backgroundColor * 255);
    alphaMap = alpha == 0;
    for i = 1:3
        tmp = squeeze(im(:,:,i));
        tmp(alphaMap) = rgbUint(i);
        im(:,:,i) = tmp;
    end
    alpha(:) = uint8(0);
end

%Resize
if (isscalar(scale) && scale ~=1) || ~isequal(scale,size(alpha))
    if matRad_cfg.isOctave
        pkg load image;
    end
    im = imresize(im,scale);
    alpha = imresize(alpha,scale);
end

end

