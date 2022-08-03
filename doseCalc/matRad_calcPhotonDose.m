function dij = matRad_calcPhotonDose(ct,stf,pln,cst,calcDoseDirect)
% matRad photon dose calculation wrapper
% 
% call
%   dij = matRad_calcPhotonDose(ct,stf,pln,cst,calcDoseDirect)
%
% input
%   ct:             ct cube
%   stf:            matRad steering information struct
%   pln:            matRad plan meta information struct
%   cst:            matRad cst struct
%   calcDoseDirect: boolian switch to bypass dose influence matrix
%                   computation and directly calculate dose; only makes
%                   sense in combination with matRad_calcDoseDirect.m
%
% output
%   dij:            matRad dij struct
%
% References
%   [1] http://www.ncbi.nlm.nih.gov/pubmed/8497215
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
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
         
    % create new engine if no engine is defined inside the pln struct
    if ~(isfield(pln, 'propDoseCalc') && isfield(pln.propDoseCalc, 'engine'))
        pln.propDoseCalc.engine = DoseEngines.matRad_DoseEnginePhotonSVD(ct,stf,pln,cst);
    end
    % if additional args are given, configure the engine
    if exist('calcDoseDirect','var')    
        pln.propDoseCalc.engine.calcDoseDirect = calcDoseDirect;
    end
    matRad_cfg.dispInfo('Starting dose calculation using %s engine.\n', pln.propDoseCalc.engine.name);
    % call calcDose from engine
    dij = matRad_calcDose(ct,stf,pln,cst);
    
end
