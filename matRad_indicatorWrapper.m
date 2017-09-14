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


% find indices of corresponding dose cubes
fieldNames = fieldnames(resultGUI);
ixCalc = []; 
% find corresponding scenario names
for k = 1:length(fieldNames)  
   res = regexp(fieldNames{k,1},[pln.bioParam.quantityVis]);
   if ~isempty(res)
      ixCalc = [ixCalc; k];
   end
end
        
for Scen = 1:numel(ixCalc)

    doseCube = resultGUI.(fieldNames{ixCalc(Scen),1});
    dvh = matRad_calcDVH(cst,doseCube,'cum');
    qi = matRad_calcQualityIndicators(cst,pln,doseCube,refGy,refVol,param);
    
    for i = 1:size(cst,1)
        cst{i,9}{Scen} = dvh{i};
        cst{i,10}{Scen} = qi{i};
    end

end

end % eof

