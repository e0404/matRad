function cst = matRad_generateBodyContour(ct,cst,thresholdHU)
% function to create a BODY contour for imported patient cases that do not have one
% 
% call
%   cst = matRad_generateBodyContour(ct,cst,thresholdHU)
%
% input
%   ct:             matrad ct structure
%   cst:            matrad cst structure
%   thresholdHU:    HU thresholding value (optional) default = -500 HU
%                  
%
% output
%   cst:        matRads cst struct  with inserted BODY contour
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();
% dicom import needs image processing toolbox -> check if available
available = matRad_checkEnvDicomRequirements(matRad_cfg.env);

if ~available
    matRad_cfg.dispError('Image processing toolbox / packages not available!');
end

if nargin<3
    thresholdHU = -500;
end
% visualize Histogram
% figure,
% histogram(ct.cubeHU{1});
% hold on
% scatter(-500,100,'x')
% hold off
% xlabel('HU')
% ylabel('#')

% Thresholding on HU
mask = zeros(ct.cubeDim);
mask(ct.cubeHU{1}>thresholdHU) = 1;
list = find(mask);
filledIm= imfill(mask);

% Write to cst
pos = size(cst,1);
cst{pos+1,1}    = pos;
cst{pos+1,2}    = 'BODY';
cst{pos+1,3}    = 'OAR';
cst{pos+1,4}{1} = find(filledIm);
cst{pos+1,5}.Priority    = 99;
cst{pos+1,5}.alphaX     = 0.1;
cst{pos+1,5}.betaX      = 0.05;
cst{pos+1,5}.Visible    = 1;
cst{pos+1,5}.visibleColor = [0.7,0.3,0.1];
cst{pos+1,6}    = [];

end