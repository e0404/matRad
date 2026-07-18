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
if nargin < 5 || isempty(visMode)
    visMode = false;
else
    sequencer.visMode = visMode;
end
if nargin < 4 || isempty(dij)
    dij = [];
end

if strcmp(pln.radiationMode, 'photons')
    % Photon (MLC leaf) sequencing goes through the class-based sequencers.

    if ~isfield(pln, 'propSeq')
        pln.propSeq = struct();
    end

    if ~isfield(pln.propSeq, 'sequencer')
        pln.propSeq.sequencer = 'siochi'; % default: siochi sequencing algorithm
        matRad_cfg.dispWarning ('pln.propSeq.sequencer not specified. Using siochi leaf sequencing (default).');
    end

    if ~any(isfield(pln.propSeq, {'numLevels', 'sequencingLevel'}))
        pln.propSeq.numLevels = 5;
        matRad_cfg.dispWarning ('pln.propSeq.sequencingLevel not specified. Using 5 sequencing levels (default).');
    elseif isfield(pln.propSeq, 'sequencingLevel')
        matRad_cfg.dispDeprecationWarning('The pln.propSeq.sequencingLevel property is deprecated. Use pln.propSeq.numLevels instead!');
        pln.propSeq.numLevels = pln.propSeq.sequencingLevel;
    end

    preconditioner = matRad_getPlnOptField(pln, 'preconditioner', false);
    dynamic = matRad_getPlnOptField(pln, 'runVMAT', false);
    continuousAperture = matRad_getPlnOptField(pln, 'continuousAperture', false);

    resultGUI = matRad_sequencePhotonsClassBased(sequencer, resultGUI, stf, dij, pln, ...
                                                 dynamic, continuousAperture, preconditioner);
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

function val = matRad_getPlnOptField(pln, fieldName, defaultVal)
% Reads pln.propOpt.(fieldName), falling back to defaultVal if pln has no
% propOpt struct or the field is unset.
if isfield(pln, 'propOpt') && isfield(pln.propOpt, fieldName)
    val = pln.propOpt.(fieldName);
else
    val = defaultVal;
end
end

function resultGUI = matRad_sequencePhotonsClassBased(sequencer, resultGUI, stf, dij, pln, ...
                                                      dynamic, continuousAperture, preconditioner)
% Siochi (static or VMAT) and static Xia/Engel sequencing through the
% class-based sequencers.

matRad_cfg = MatRad_Config.instance();

% bridge legacy pln.propOpt/propSeq fields onto the sequencer object,
% since assignPropertiesFromPln only auto-maps pln.propSeq.*
sequencer.sequencingLevel = pln.propSeq.numLevels;
sequencer.preconditioner = preconditioner;
if isa(sequencer, 'matRad_PhotonSequencerVMATAbstract')
    sequencer.runVMAT = dynamic;
    sequencer.continuousAperture = continuousAperture;
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
end
