function resultGUI = matRad_sequencing(resultGUI,stf,dij,pln,visBool)
% matRad inverse planning wrapper function
% 
% call
%   resultGUI = matRad_sequencing(resultGUI,stf,dij,pln)
%
% input
%   dij:        matRad dij struct
%   stf:        matRad stf struct
%   pln:        matRad pln struct
%   resultGUI:  struct containing optimized fluence vector, dose, and (for
%               biological optimization) RBE-weighted dose etc.
%
% output
%   resultGUI:  struct containing optimized fluence vector, dose, and (for
%               biological optimization) RBE-weighted dose etc.
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team. 
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

if nargin < 5
    visBool = 0;
end

if strcmp(pln.radiationMode,'photons') && (pln.propSeq.runSequencing || pln.propOpt.runDAO)
    
    if ~isfield(pln.propSeq, 'sequencer')
        pln.propSeq.sequencer = 'siochi'; % default: siochi sequencing algorithm
        matRad_cfg.dispWarning ('pln.propSeq.sequencer not specified. Using siochi leaf sequencing (default).')
    end
    
    if ~isfield(pln.propSeq, 'sequencingLevel')
        pln.propSeq.sequencingLevel = 5;
         matRad_cfg.dispWarning ('pln.propSeq.sequencingLevel not specified. Using 5 sequencing levels (default).')
    end
    
    switch pln.propSeq.sequencer
        case 'xia'
            resultGUI = matRad_xiaLeafSequencing(resultGUI,stf,dij,pln.propSeq.sequencingLevel,visBool);
        case 'engel'
            resultGUI = matRad_engelLeafSequencing(resultGUI,stf,dij,pln.propSeq.sequencingLevel,visBool);
        case 'siochi'
            resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln.propSeq.sequencingLevel,visBool);
        otherwise
            matRad_cfg.dispError('Could not find specified sequencing algorithm');
    end
elseif (pln.propSeq.runSequencing || pln.propOpt.runDAO) && ~strcmp(pln.radiationMode,'photons')
    matRad_cfg.dispWarning('Sequencing is only specified for pln.radiationMode = "photons". Continuing with out sequencing ... ')
end
end


