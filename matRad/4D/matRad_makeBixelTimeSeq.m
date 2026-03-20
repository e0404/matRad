function timeSequence = matRad_makeBixelTimeSeq(stf, resultGUI)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% using the steering information of matRad, makes a time sequenced order
% according to the irradiation scheme in spot scanning
%
% call
%   timeSequence = matRad_makeBixelTimeSeq(stf, resultGUI)
%
% input
%   stf:            matRad steering information struct
%   resultGUI:      struct containing optimized fluence vector
%
% output
%   timeSequence:      struct containing bixel ordering information and the
%                   time sequence of the spot scanning
%
% References
%   spill structure and timing informations:
%   http://cdsweb.cern.ch/record/1182954
%   http://iopscience.iop.org/article/10.1088/0031-9155/56/20/003/meta
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();
matRad_cfg.dispWarning('This function is Outdated use new SequencingClass');

sequencer = matRad_ParticleSequencer();

timeSequence = sequencer.sequence(stf, resultGUI.w);

end
