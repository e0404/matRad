function dij = matRad_calcPhotonDoseMC(ct,stf,pln,cst,nCasePerBixel,visBool)
% matRad ompMC monte carlo photon dose calculation wrapper
%
% call
%   dij = matRad_calcPhotonDoseMc(ct,stf,pln,cst,visBool)
%
% input
%   ct:                         matRad ct struct
%   stf:                        matRad steering information struct
%   pln:                        matRad plan meta information struct
%   cst:                        matRad cst struct
%   visBool:                    binary switch to enable visualization
% output
%   dij:                        matRad dij struct
%
% References
%   -
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


matRad_cfg =  MatRad_Config.instance();

matRad_cfg.dispDeprecationWarning('The old dose calculation functions are deprecated! Try to use matRad_calcDoseInfluence with the new engine format from now on!');

% could be also set as pln property e.g pln.propDoseCalc.useDeprecated
if isfield(pln, 'propDoseCalc') && isfield(pln.propDoseCalc, 'engine')
    matRad_cfg.dispWarning('You should not use the deprecated MC calculation with the new engine architecture! Setting ompMC as engine!');
end

engine = DoseEngines.matRad_PhotonOmpMCEngine(pln);

% assign old deprecated defaults
if exist('nCasePerBixel','var')
    engine.numHistoriesPerBeamlet = nCasePerBixel;
else
    engine.numHistoriesPerBeamlet = matRad_cfg.propMC.ompMC_defaultHistories;
end
engine.outputMCvariance = matRad_cfg.propMC.ompMC_defaultOutputVariance;
    
if exist('visBool','var')
    engine.visBool = visBool;
end

matRad_cfg.dispInfo('Starting dose calculation using %s engine.\n', engine.name);

pln.propDoseCalc = engine;

% call the calcDose from engine
dij = matRad_calcDoseInfluence(ct,cst,stf,pln);
end