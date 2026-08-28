function fluenceObjectives = matRad_getFluenceObjectives(pln, cst, dij)
% matRad helper collecting and initializing the fluence objectives of a plan
%
% Fluence objectives (see matRad_FluenceOptimizationFunction) act on the
% bixel weight vector directly instead of on the dose. They can be given in
% two places::
%
%   pln.propOpt.fluenceObjectives   cell array (or single object) of
%                                   objectives that are not tied to a
%                                   particular structure - the usual place
%   cst{i,6}                        alongside the dose objectives of a
%                                   structure; the structure itself is
%                                   ignored, the objective still acts on the
%                                   whole fluence
%
% A fluence objective is a pure specification until it is set up for a
% beamlet geometry. This function performs that setup from the dij, so that
% the caller only has to hand over the collected objectives.
%
% call:
%   fluenceObjectives = matRad_getFluenceObjectives(pln,cst,dij)
%
% input:
%   pln:    matRad pln struct
%   cst:    matRad cst struct
%   dij:    matRad dij struct, providing the beamlet geometry the objectives
%           are set up for. Pass [] to only collect them without setup.
%
% output:
%   fluenceObjectives:  cell array of initialized fluence objectives
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2026 the matRad development team.
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

fluenceObjectives = {};

% from pln.propOpt
if isfield(pln, 'propOpt') && isfield(pln.propOpt, 'fluenceObjectives')
    plnObjectives = pln.propOpt.fluenceObjectives;
    if ~iscell(plnObjectives)
        plnObjectives = {plnObjectives};
    end
    fluenceObjectives = [fluenceObjectives, plnObjectives(:)'];
end

% from the cst - dose and fluence functions can be mixed there, the dose
% objective loops filter by class anyway
for i = 1:size(cst, 1)
    for j = 1:numel(cst{i, 6})
        if isa(cst{i, 6}{j}, 'matRad_FluenceOptimizationFunction')
            fluenceObjectives{end + 1} = cst{i, 6}{j}; %#ok<AGROW>
        end
    end
end

if isempty(fluenceObjectives) || isempty(dij)
    return
end

% set up for this dij's beamlet geometry and validate
for i = 1:numel(fluenceObjectives)
    obj = fluenceObjectives{i};

    if ~isa(obj, 'matRad_FluenceOptimizationFunction')
        matRad_cfg.dispError(['Entry %d of the fluence objectives is not a ' ...
                              'matRad_FluenceOptimizationFunction!'], i);
    end

    obj = obj.setupForDij(dij);

    [isCompat, msg] = obj.isCompatible(dij.totalNumOfBixels);
    if ~isCompat
        matRad_cfg.dispError('%s', msg);
    end

    fluenceObjectives{i} = obj;
end

matRad_cfg.dispInfo('Using %d fluence objective(s) in addition to the dose objectives\n', numel(fluenceObjectives));
