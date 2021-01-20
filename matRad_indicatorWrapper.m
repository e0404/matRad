function [dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUI,refGy,refVol)
% matRad indictor wrapper
% 
% call
%   [dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUI)
%   [dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUI,refGy,refVol)
%
% input
%   cst:                  matRad cst struct
%   pln:                  matRad pln struct
%   resultGUI:            matRad resultGUI struct
%   refGy: (optional)     array of dose values used for V_XGy calculation
%                         default is [40 50 60]
%   refVol:(optional)     array of volumes (0-100) used for D_X calculation
%                         default is [2 5 95 98]
%                         NOTE: Call either both or none!
%
% output
%   dvh: matRad dvh result struct
%   qi:  matRad quality indicator result struct
%   graphical display of all results
%
% References
%   -
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

if isfield(resultGUI,'RBExDose')
    doseCube = resultGUI.RBExDose;
else
    doseCube = resultGUI.physicalDose;
end

if ~exist('refVol', 'var') 
    refVol = [];
end

if ~exist('refGy', 'var')
    refGy = [];
end

dvh = matRad_calcDVH(cst,doseCube,'cum');
qi  = matRad_calcQualityIndicators(cst,pln,doseCube,refGy,refVol);

figure,set(gcf,'Color',[1 1 1]);
subplot(2,1,1)
matRad_showDVH(dvh,cst,pln);
subplot(2,1,2)
ixVoi = cellfun(@(c) c.Visible == 1,cst(:,5));
qi = qi(ixVoi);
matRad_showQualityIndicators(qi);



