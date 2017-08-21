function cst = matRad_indicatorWrapper(cst,pln,resultGUI,refGy,refVol,param)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad indictor wrapper
% 
% call
%   matRad_calcIndicators(cst,pln,cube,dvhType,param,refGy,refVol,lineStyleIndicator)
%
% input
%   cst:                  matRad cst struct
%   pln:                  matRad pln struct
%   resultGUI:            matRad resultGUI struct
%
% output
%   various quality indicators as well as dvh stored in cst
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team. 
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

if exist('param','var')
   if ~isfield(param,'logLevel')
      param.logLevel = 1;
   end
else
   param.logLevel = 1;
end


dvh = matRad_calcDVH(cst,doseCube,'cum');
qi = matRad_calcQualityIndicators(cst,pln,doseCube,refGy,refVol,param);

numOfScenarios = 1;
for i = 1:size(cst,1)
    % overload with scenarios
    cst{i,8} = cell(numOfScenarios,1);
    cst{i,9} = cell(numOfScenarios,1);
    
    cst{i,8}{1} = dvh{i};
    cst{i,9}{1} = qi{i};
end    


end % eof

