classdef (Abstract) matRad_DoseEngineMonteCarlo < DoseEngines.matRad_DoseEngine
% matRad_DoseEngineMonteCarlo: abstract superclass for all dose calculation 
%   engines which are based on monte carlo calculation 
%   for more informations see superclass
%   DoseEngines.matRad_DoseEngine
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (SetAccess = public, GetAccess = public)
                         
        nCasePerBixel ; %number of histories per beamlet
        
    end
        
    methods
        
        function this = matRad_DoseEngineMonteCarlo()
            
            % call superclass constructor
            this = this@DoseEngines.matRad_DoseEngine();
   
            % create config instance
            matRad_cfg = MatRad_Config.instance();
            
            %set number of particles simulated per pencil beam
            this.nCasePerBixel = matRad_cfg.propMC.MCsquare_defaultHistories;
            
        end
    end
    
end

