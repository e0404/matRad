function resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst,w,nHistories)
    % matRad function to bypass dij calculation
    %   Should not be used directly anymore as it is deprecated. Use matRad_calcDoseForward instead.
    % 
    % call
    %   resultGUI = matRad_calcDoseDirec(ct,stf,pln,cst)
    %   resultGUI = matRad_calcDoseDirec(ct,stf,pln,cst,w)
    %
    % input
    %   ct:         ct cube
    %   stf:        matRad steering information struct
    %   pln:        matRad plan meta information struct
    %   cst:        matRad cst struct
    %   w:          (optional, if no weights available in stf): bixel weight
    %               vector
    %
    % output
    %   resultGUI:  matRad result struct
    %
    % References
    %   -
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2024 the matRad development team. 
    % 
    % This file is part of the matRad project. It is subject to the license 
    % terms in the LICENSE file found in the top-level directory of this 
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
    % of the matRad project, including this file, may be copied, modified, 
    % propagated, or distributed except according to the terms contained in the 
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    matRad_cfg =  MatRad_Config.instance();
    matRad_cfg.dispDeprecationWarning('This function is deprecated. Please use matRad_calcDoseForward with appropriate engine set in pln.propDoseCalc.');
    
    % dose calculation
    switch pln.radiationMode
        case {'protons','helium','carbon'}
            engine = DoseEngines.matRad_ParticleHongPencilBeamEngine(pln);
        case 'photons'
            engine = DoseEngines.matRad_PhotonPencilBeamSVDEngine(pln);
        otherwise
            matRad_cfg.dispError('Radiation mode ''%s'' not supported!',pln.radiationMode)
    end
    
    if nargin < 4
        resultGUI = engine.calcDoseForward(ct,cst,stf);
    else
        resultGUI = engine.calcDoseForward(ct,cst,stf,w);
    end
    
    
    
    
    
    
    