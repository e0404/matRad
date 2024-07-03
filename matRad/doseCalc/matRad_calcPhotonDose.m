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

matRad_cfg = MatRad_Config.instance();
matRad_cfg.dispDeprecationWarning('The old dose calculation functions are deprecated! Try to use matRad_calcDoseInfluence with the new engine format from now on!');

% could be also set as pln property e.g pln.propDoseCalc.useDeprecated
if isfield(pln, 'propDoseCalc') && isfield(pln.propDoseCalc, 'engine')
    matRad_cfg.dispWarning('You should not use the deprecated MC calculation with the new engine architecture! Setting SVD pencil beam as engine!');
end

engine = DoseEngines.matRad_PhotonPencilBeamSVDEngine(pln);

% if additional args are given, configure the engine
if exist('calcDoseDirect','var')
    engine.calcDoseDirect = calcDoseDirect;
end
matRad_cfg.dispInfo('Starting dose calculation using %s engine.\n', engine.name);

pln.propDoseCalc = engine;
% call calcDose from engine
dij = matRad_calcDoseInfluence(ct,cst,stf,pln);
    
end
