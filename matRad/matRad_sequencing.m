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
    matRad_cfg.dispDeprecationWarning(['The argument order of matRad_sequencing changed to ' ...
                                       'matRad_sequencing(resultGUI, stf, pln, dij). Please update your call.']);
    [pln, dij] = deal(dij, pln);
end

sequencer = matRad_SequencerBase.getSequencerFromPln(pln);

% Handle optional inputs
if nargin >= 5 && ~isempty(visMode)
    sequencer.visMode = visMode;
end
if nargin < 4 || isempty(dij)
    dij = [];
end

if strcmp(pln.radiationMode, 'photons')
    % Photon (MLC leaf) sequencing goes through the class-based sequencers.
    % The canonical sequencer configuration lives under pln.propSeq and was
    % already auto-mapped onto the sequencer during its construction; here
    % only the deprecated pln.propOpt locations are bridged (honored for one
    % release) plus the settings that genuinely live outside propSeq.

    if ~isfield(pln, 'propSeq')
        pln.propSeq = struct();
    end
    if ~isfield(pln, 'propOpt')
        pln.propOpt = struct();
    end

    if ~isfield(pln.propSeq, 'preconditioner') && isfield(pln.propOpt, 'preconditioner')
        matRad_cfg.dispDeprecationWarning('pln.propOpt.preconditioner is deprecated. Use pln.propSeq.preconditioner instead!');
        sequencer.preconditioner = pln.propOpt.preconditioner;
    end

    % pln.propOpt.runVMAT deliberately remains the canonical mode flag
    % (like runDAO): optimization needs it first to select the problem
    % class, so it is bridged here rather than moved to propSeq.
    dynamic = matRad_getFieldOrDefault(pln.propOpt, 'runVMAT', false);

    if isa(sequencer, 'matRad_PhotonSequencerVMATAbstract')
        sequencer.runVMAT = dynamic;
        if ~isfield(pln.propSeq, 'continuousAperture') && isfield(pln.propOpt, 'continuousAperture')
            matRad_cfg.dispDeprecationWarning('pln.propOpt.continuousAperture is deprecated. Use pln.propSeq.continuousAperture instead!');
            sequencer.continuousAperture = pln.propOpt.continuousAperture;
        end
        if ~isempty(dij) && isfield(dij, 'weightToMU')
            sequencer.weightToMU = dij.weightToMU;
        end
    elseif dynamic
        matRad_cfg.dispError(['Sequencer ''%s'' does not support VMAT (dynamic) delivery. ' ...
                              'Use ''siochi'' for VMAT plans.'], sequencer.shortName);
    end

    sequence = sequencer.sequence(resultGUI.w, stf);
    if ~isempty(dij)
        % merge the computed dose cubes into resultGUI instead of overwriting
        % it, so that pre-existing fields (e.g. from FMO) survive
        doseCubes = matRad_calcCubes(sequence.w, dij);
        fNames = fieldnames(doseCubes);
        for f = 1:numel(fNames)
            resultGUI.(fNames{f}) = doseCubes.(fNames{f});
        end
    else
        matRad_cfg.dispWarning('Dose not recalculated with sequenced fluence');
    end
    resultGUI.sequencing = sequence;
    if isfield(sequence, 'apertureInfo')
        resultGUI.apertureInfo = sequence.apertureInfo;
    end
else
    % Non-photon (e.g. particle) sequencing goes through the new
    % class-based sequencer, which only derives spot delivery order/timing
    % and leaves the fluence - and thus the dose - unchanged.
    sequence = sequencer.sequence(resultGUI.w, stf);
    resultGUI.sequencing = sequence;
    if isfield(sequence, 'apertureInfo')
        resultGUI.apertureInfo = sequence.apertureInfo;
    end
end

end
