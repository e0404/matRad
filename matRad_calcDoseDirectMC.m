function resultGUI = matRad_calcDoseDirectMC(ct,stf,pln,cst,w,nHistories)
% matRad function to bypass dij calculation for MC dose calculation 
% matRad dose calculation wrapper for MC dose calculation algorithms
% bypassing dij calculation for MC dose calculation algorithms.
% 
% call
%   resultGUI = matRad_calcDoseDirecMC(ct,stf,pln,cst)
%   resultGUI = matRad_calcDoseDirecMC(ct,stf,pln,cst,w)
%   resultGUI = matRad_calcDoseDirectMC(ct,stf,pln,cst,nHistories)
%   resultGUI = matRad_calcDoseDirectMC(ct,stf,pln,cst,w,nHistories)
%
% input
%   ct:         ct cube
%   stf:        matRad steering information struct
%   pln:        matRad plan meta information struct
%   cst:        matRad cst struct
%   w:          (optional, if no weights available in stf): bixel weight
%               vector
%   nHistories: (optional) number of histories
%
% output
%   resultGUI:  matRad result struct
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


matRad_cfg =  MatRad_Config.instance();
matRad_cfg.dispDeprecationWarning('This function is deprecated. Please use matRad_calcDoeDirect with appropriate Monte Carlo engine set in pln.propDoseCalc.');

% dose calculation
switch pln.radiationMode
    case 'protons'
        engine = DoseEngines.matRad_ParticleMCsquareEngine(pln);
    case 'photons'
        engine = DoseEngines.matRad_PhotonOmpMCEngine(pln);
    otherwise
        matRad_cfg.dispError('Radiation mode ''%s'' not supported!',pln.radiationMode)
end

if nargin == 5
    engine.numHistoriesDirect = nHistories;
end

if nargin < 4
    resultGUI = engine.calcDoseForward(ct,cst,stf);
else
    resultGUI = engine.calcDoseForward(ct,cst,stf,w);
end






