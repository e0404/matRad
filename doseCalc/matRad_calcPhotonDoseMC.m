function dij = matRad_calcPhotonDoseMC(ct,stf,pln,cst,nCasePerBixel,visBool,useDeprecated)
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
%   useDeprecated:              optional var for using the "old" function, when false
%                               use dose engines
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

    % could be also set as pln property e.g pln.propDoseCalc.useDeprecated
    if ~(isfield(pln, 'propDoseCalc') && isfield(pln.propDoseCalc, 'engine'))
        % create new engine if no engine is defined inside the pln struct
        pln.propDoseCalc.engine = DoseEngines.matRad_PhotonMonteCarloEngineOmpMC(ct,stf,pln,cst);
    end
    % if additional args are given, configure the engine
    if exist('nCasePerBixel','var')
        pln.propDoseCalc.engine.nCasePerBixel = nCasePerBixel;
    end
    if exist('visBool','var')    
        pln.propDoseCalc.engine.visBool = visBool;
    end
    matRad_cfg.dispInfo('Starting dose calculation using %s engine.\n', pln.propDoseCalc.engine.name);
    
    % call the calcDose from engine
    dij = matRad_calcDose(ct,stf,pln,cst);
end