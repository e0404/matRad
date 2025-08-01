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
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sequencer = matRad_SequencingBase.getSequencerFromPln(pln);

if nargin == 5
    sequencer.visMode = visMode;

end

sequence = sequencer.sequence(resultGUI.w,stf);
resultGUI = sequencer.updateResultGUI(resultGUI,sequence,stf,dij);

end


