function dij = matRad_calcParticleDoseMC(ct,stf,pln,cst,nHistories,calcDoseDirect)
% matRad Monte Carlo particle dose calculation wrapper
%   Will call the appropriate subfunction for the respective
%   MC dose-calculation engine
%
% call
%   dij = matRad_calcParticleDoseMc(ct,stf,pln,cst)
%   dij = matRad_calcParticleDoseMc(ct,stf,pln,cst,calcDoseDirect)
%
% input
%   ct:                         matRad ct struct
%   stf:                        matRad steering information struct
%   pln:                        matRad plan meta information struct
%   cst:                        matRad cst struct
%   nHistories:                 number of histories per beamlet
%   calcDoseDirect:             (optional) binary switch to enable forward 
%                               dose calcualtion (default: false)
% output
%   dij:                        matRad dij struct
%
% References
%
% -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team. 
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

if nargin < 6
    calcDoseDirect = false;
end

if ~isfield(pln,'propMC') || ~isfield(pln.propMC,'proton_engine')
    matRad_cfg.dispInfo('Using default proton MC engine "%s"\n',matRad_cfg.propMC.default_proton_engine);
    engine = matRad_cfg.propMC.default_proton_engine;
else
    engine = pln.propMC.proton_engine;
end


if nargin < 5 
    if ~calcDoseDirect
        nHistories = matRad_cfg.propMC.particles_defaultHistories;
        matRad_cfg.dispInfo('Using default number of Histories per bixel: %d\n',nHistories);
    else
        nHistories = matRad_cfg.propMC.direct_defaultHistories;
        matRad_cfg.dispInfo('Using default number of Histories for forward dose calculation: %d\n',nHistories);
    end
end

switch engine
    case 'MCsquare'
        dij = matRad_calcParticleDoseMCsquare(ct,stf,pln,cst,nHistories,calcDoseDirect);
    case 'TOPAS'
        dij = matRad_calcParticleDoseMCtopas(ct,stf,pln,cst,nHistories,calcDoseDirect);
    otherwise
        matRad_cfg.dispError('MC engine %s not known/supported!',engine);
end

end
