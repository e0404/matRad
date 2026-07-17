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
if nargin < 5 || isempty(visMode)
    visMode = false;
else
    sequencer.visMode = visMode;
end
if nargin < 4 || isempty(dij)
    dij = [];
end

if strcmp(pln.radiationMode, 'photons')
    % Photon (MLC leaf) sequencing keeps using the functional
    % implementations, which carry VMAT-specific dynamic/arc-sequencing
    % support that the class-based sequencers do not implement yet.

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

    if isfield(pln, 'propOpt') && isfield(pln.propOpt, 'preconditioner')
        preconditioner = pln.propOpt.preconditioner;
    else
        preconditioner = false;
    end

    if isfield(pln, 'propOpt') && isfield(pln.propOpt, 'runVMAT')
        dynamic = pln.propOpt.runVMAT;
    else
        dynamic = false;
    end

    if isfield(pln, 'propOpt') && isfield(pln.propOpt, 'numApertures')
        numApertures = pln.propOpt.numApertures;
    else
        numApertures = 0;
    end

    if isfield(pln, 'propOpt') && isfield(pln.propOpt, 'continuousAperture')
        continuousAperture = pln.propOpt.continuousAperture;
    else
        continuousAperture = false;
    end

    varArgList = { ...
                  'visBool', visMode, ...
                  'dynamic', dynamic, ...
                  'numApertures', numApertures, ...
                  'continuousAperture', continuousAperture, ...
                  'preconditioner', preconditioner};

    % Could probably consolidate a lot of the code in the following
    % functions.
    switch pln.propSeq.sequencer
        case 'xia'
            resultGUI = matRad_xiaLeafSequencing(resultGUI, stf, dij, pln.propSeq.numLevels, varArgList{:});
        case 'engel'
            resultGUI = matRad_engelLeafSequencing(resultGUI, stf, dij, pln.propSeq.numLevels, varArgList{:});
        case 'siochi'
            resultGUI = matRad_siochiLeafSequencing(resultGUI, stf, dij, pln.propSeq.numLevels, varArgList{:});
        otherwise
            matRad_cfg.dispError('Could not find specified sequencing algorithm ''%s''', pln.propSeq.sequencer);
    end

    % keep the aperture info available under resultGUI.sequencing.apertureInfo
    % too, matching the location used by the class-based sequencers
    if isfield(resultGUI, 'apertureInfo') && isfield(resultGUI, 'sequencing') && isstruct(resultGUI.sequencing)
        resultGUI.sequencing.apertureInfo = resultGUI.apertureInfo;
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
