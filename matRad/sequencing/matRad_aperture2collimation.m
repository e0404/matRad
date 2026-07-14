function [pln, stf] = matRad_aperture2collimation(pln, stf, sequencing, apertureInfo)
% matRad function to convert sequencing information / aperture information
% into collimation information in pln and stf for field-based dose
% calculation
% 
% call:
%   [pln,stf] = matRad_aperture2collimation(pln,stf,sequencing,apertureInfo)
%
% input:
%   pln:            pln file used to generate the sequenced plan
%   stf:            stf file used to generate the sequenced plan
%   sequencing:     sequencing information (from resultGUI)
%   apertureInfo:   apertureInfo (from resultGUI)
%
% output:
%   pln:            matRad pln struct with collimation information
%   stf:            matRad stf struct with shapes instead of beamlets
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022-2026 the matRad development team.
% Author: wahln
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

matRad_cfg.dispWarning('This function is outdated, use the sequencer class instead');

sequencer = matRad_SequencerBase.getSequencerFromPln(pln);

if ~isfield(sequencing, 'apertureInfo')
    if nargin < 4
        matRad_cfg.dispError('Sequencing struct does not contain apertureInfo and none was provided!');
    end
    sequencing.apertureInfo = apertureInfo;
end
[pln, stf] = sequencer.aperture2collimation(pln, stf, sequencing);

end
