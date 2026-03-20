function timeSequence = matRad_makePhaseMatrix(timeSequence, numOfPhases, motionPeriod, motion)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   using the time sequence and the ordering of the bixel iradiation, and
%   number of scenarios, makes a phase matrix of size number of bixels *
%   number of scenarios
%
%
% call
%   timeSequence = matRad_makePhaseMatrix(timeSequence, numOfPhases, motionPeriod, motion)
%
% input
%   timeSequence:   struct containing bixel ordering information and the
%                   time sequence of the spot scanning
%   numOfCtScen:    number of the desired phases
%   motionPeriod:   the extent of a whole breathing cycle (in seconds)
%   motion:         motion scenario: 'linear'(default), 'sampled_period'
%
% output
%   timeSequence:      phase matrix field added
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
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();
matRad_cfg.dispWarning('This function is Outdated use new SequencingClass');

sequencer = matRad_ParticleSequencer();
timeSequence = sequencer.makePhaseMatrix(timeSequence, numOfPhases, motionPeriod);

end
