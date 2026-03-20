function resultGUI = matRad_engelLeafSequencing(resultGUI, stf, dij, numOfLevels, visBool)
% multileaf collimator leaf sequencing algorithm
% for intensity modulated beams with multiple static segments accroding
% to Engel et al. 2005 Discrete Applied Mathematics
%
% call
%   resultGUI = matRad_engelLeafSequencing(resultGUI,stf,dij,numOfLevels,visBool)
%
% input
%   resultGUI:          resultGUI struct to which the output data will be added, if
%                       this field is empty resultGUI struct will be created
%   stf:                matRad steering information struct
%   dij:                matRad's dij matrix
%   numOfLevels:        number of stratification levels
%   visBool:            toggle on/off visualization (optional)
%
% output
%   resultGUI:          matRad result struct containing the new dose cube
%                       as well as the corresponding weights
%
% References
%   [1] http://www.sciencedirect.com/science/article/pii/S0166218X05001411
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matRad_cfg = MatRad_Config.instance();

matRad_cfg.dispWarning('This function is Outdated use new SequencingClass');

sequencer = matRad_SequencingPhotonsEngelLeaf();
sequencer.visMode = visBool;
sequencer.sequencingLevel = numOfLevels;

sequence = sequencer.sequence(resultGUI.w, stf);
if ~isempty(dij)
    resultGUI = matRad_calcCubes(sequence.w, dij);
else
    matRad_cfg.dispWarning('Dose not recalcaulted with sequenced fluence');
end
resultGUI.sequencing   = sequence;

end
