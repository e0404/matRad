function dij = matRad_calcDoseInfluence(ct,cst,stf,pln)
% matRad dose calculation automaticly creating the appropriate dose engine
% for the given pln struct and called the associated dose calculation funtion
%
% call
%   dij =  matRad_calcDoseInfluence(ct,cst,stf,pln)
%
% input
%   ct:         ct cube
%   cst:        matRad cst cell array
%   stf:        matRad steering information struct
%   pln:        matRad plan meta information struct
%
%
% output
%   dij:            matRad dij struct
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
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


matRad_cfg = MatRad_Config.instance();

%Deprecation warnings
if ~isfield(stf,'machine')
    matRad_cfg.dispDeprecationWarning('stf should contain the machine name in the ''machine'' field since matRad 3. Manually adding ''%s'' from pln.',pln.machine);
    for i=1:numel(stf)
        stf(i).machine = pln.machine;
    end
end

engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);

%call the calcDose funktion
dij = engine.calcDoseInfluence(ct,cst,stf);

end
