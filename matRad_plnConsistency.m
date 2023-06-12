function pln = matRad_plnConsistency(pln)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_plnConsistency checks the consistency of the optimized quantities 
% for the multiple modalities; all optimization quantities should be the
% same!
%
% call
%   pln = matRad_plnConsistency(pln, param)
%
% input
%   pln(:):          pln stuct for multiple modalities  
% output
%   pln struct
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

matRad_cfg = MatRad_Config.instance();
if numel(pln) == 1
    matRad_cfg.dispInfo(['single pln is by definition consistent with itself... \n']);
elseif numel(pln) > 1
    maxh = [];
    % find level in the heirarchy (physical, RBExD, Effect)
    for i = 1:numel(pln)
        switch pln(i).bioParam.quantityOpt
            case{'physicalDose'}
                maxh = [maxh 1];
            case{'RBExD'}
                maxh = [maxh 2];
            case{'effect'}
                maxh = [maxh 3];
        end
    end
    
    if numel(unique(maxh)) > 1
        maxh = max(maxh);
        %This will probably clash with new class based bioModel implementation
        quantityOpt = pln(i).bioParam.AvailableQuantitiesForOpt{maxh};
        matRad_cfg.dispInfo(['make plan consistent. optimized quantity chosen by hirarchy is ' quantityOpt '. \n' ]);
        for i = 1:numel(pln)
            pln(i).bioParam = matRad_bioModel(pln(i).radiationMode, quantityOpt, pln(i).bioParam.model);
        end
    elseif numel(unique(maxh)) == 1
        matRad_cfg.dispInfo(['pln is consistent... \n']);
    end
end
