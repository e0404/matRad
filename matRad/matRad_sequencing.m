function resultGUI = matRad_sequencing(resultGUI, stf, pln, dij, visMode)
% matRad inverse planning wrapper function
% 
% call:
%   resultGUI = matRad_sequencing(resultGUI,stf,pln,dij)
%   resultGUI = matRad_sequencing(resultGUI,stf,pln,dij,visMode)
%
% input:
%   resultGUI:  struct containing optimized fluence vector, dose, and (for
%               biological optimization) RBE-weighted dose etc.
%   stf:        matRad stf struct
%   pln:        matRad pln struct
%   dij:        matRad dij struct (optional; if given, the dose is
%               recomputed from the sequenced fluence for photon plans)
%   visMode:    toggle sequencing visualization on/off (optional)
%
% output:
%   resultGUI:  struct containing optimized fluence vector, dose, and (for
%               biological optimization) RBE-weighted dose etc.
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016-2026 the matRad development team.
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

% Backwards compatibility: the order of the pln and dij arguments was swapped
% (previously matRad_sequencing(resultGUI, stf, dij, pln)). A dij always carries
% the dose-influence grid (doseGrid), a pln never does at the top level - so if
% the pln position holds a dij, the call uses the old order and the arguments
% are swapped with a deprecation warning.
if nargin >= 4 && isstruct(pln) && isfield(pln, 'doseGrid')
    matRad_cfg.dispDeprecationWarning('The argument order of matRad_sequencing changed to matRad_sequencing(resultGUI, stf, pln, dij). Please update your call.');
    [pln, dij] = deal(dij, pln);
end

sequencer = matRad_SequencerBase.getSequencerFromPln(pln);

% Handle optional inputs
if nargin == 5 && ~isempty(visMode)
    sequencer.visMode = visMode;
end
if nargin < 4 || isempty(dij)
    dij = [];
end

sequence = sequencer.sequence(resultGUI.w, stf);

% Aperture-based (photon) sequencing modifies the fluence into deliverable
% MLC segments, so the dose has to be recomputed from the sequenced fluence.
% Particle sequencing only derives the spot delivery order/timing, leaves the
% fluence unchanged and returns a per-beam struct array - the existing dose
% cubes stay valid and must not be recomputed here.
if isa(sequencer, 'matRad_PhotonSequencerAbstract')
    if ~isempty(dij)
        resultGUI = matRad_calcCubes(sequence.w, dij);
    else
        matRad_cfg.dispWarning('Dose not recalculated with sequenced fluence');
    end
end
resultGUI.sequencing   = sequence;

% keep a backward-compatible copy of the aperture info at the top level so
% that legacy calls (e.g. matRad_directApertureOptimization) still work
if isfield(sequence, 'apertureInfo')
    resultGUI.apertureInfo = sequence.apertureInfo;
end

end
