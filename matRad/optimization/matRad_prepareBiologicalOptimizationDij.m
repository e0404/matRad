function dij = matRad_prepareBiologicalOptimizationDij(dij, cst, backProjection)
% matRad_prepareBiologicalOptimizationDij prepares biological dij metadata
%
% call
%   dij = matRad_prepareBiologicalOptimizationDij(dij,cst,backProjection)
%
% input
%   dij:              matRad dij struct
%   cst:              matRad cst struct on the dose grid
%   backProjection:   matRad backprojection object
%
% output
%   dij:              dij with biological ax/bx, ixDose and gamma metadata
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

if ~isa(backProjection, 'matRad_EffectProjection')
    return
end

matRad_cfg = MatRad_Config.instance();
numOfVoxels = matRad_getDoseGridNumOfVoxels(dij);
numOfCtScen = matRad_getNumOfCtScenarios(dij, cst);

axBxFromCst = ~all(isfield(dij, {'ax', 'bx'}));
if axBxFromCst
    matRad_cfg.dispWarning(['Biological optimization requested, but no ' ...
                            'ax & bx provided in dij. Getting from cst...']);
    [dij.ax, dij.bx] = matRad_getPhotonLQMParameters(cst, numOfVoxels);
    ixZeroDose = matRad_getZeroDoseVoxelMask(dij, numOfVoxels);
    for scenIx = 1:numel(dij.ax)
        dij.ax{scenIx}(ixZeroDose) = 0;
        dij.bx{scenIx}(ixZeroDose) = 0;
    end
else
    dij.ax = matRad_asCtScenarioCell(dij.ax, numOfCtScen, ...
                                     numOfVoxels, 'dij.ax');
    dij.bx = matRad_asCtScenarioCell(dij.bx, numOfCtScen, ...
                                     numOfVoxels, 'dij.bx');
end

[dij.ax, dij.bx] = matRad_validateAxBxCells(dij.ax, dij.bx, ...
                                            numOfCtScen, numOfVoxels);
if ~axBxFromCst
    matRad_validateAxBxAgainstCst(dij.ax, dij.bx, cst, numOfVoxels);
end
dij.ixDose = cell(numOfCtScen, 1);
for scenIx = 1:numOfCtScen
    dij.ixDose{scenIx} = dij.bx{scenIx} ~= 0;
end

if isa(backProjection, 'matRad_VariableRBEProjection')
    dij.gamma = cell(numOfCtScen, 1);
    for scenIx = 1:numOfCtScen
        dij.gamma{scenIx} = zeros(numOfVoxels, 1);
        dij.gamma{scenIx}(dij.ixDose{scenIx}) = ...
            dij.ax{scenIx}(dij.ixDose{scenIx}) ./ ...
            (2 * dij.bx{scenIx}(dij.ixDose{scenIx}));
    end
end

end

function numOfVoxels = matRad_getDoseGridNumOfVoxels(dij)
matRad_cfg = MatRad_Config.instance();
if ~isstruct(dij) || ~isfield(dij, 'doseGrid') || ...
        ~isfield(dij.doseGrid, 'numOfVoxels') || ...
        ~isnumeric(dij.doseGrid.numOfVoxels) || ...
        ~isscalar(dij.doseGrid.numOfVoxels) || ...
        dij.doseGrid.numOfVoxels < 1
    matRad_cfg.dispError(['dij.doseGrid.numOfVoxels must be a positive ' ...
                          'numeric scalar for biological optimization!']);
end
numOfVoxels = dij.doseGrid.numOfVoxels;
end

function numOfCtScen = matRad_getNumOfCtScenarios(dij, cst)
matRad_cfg = MatRad_Config.instance();
if isfield(dij, 'physicalDose') && iscell(dij.physicalDose) && ...
        ~isempty(dij.physicalDose)
    numOfCtScen = size(dij.physicalDose, 1);
elseif ~isempty(cst) && size(cst, 2) >= 4 && iscell(cst{1, 4})
    numOfCtScen = numel(cst{1, 4});
else
    matRad_cfg.dispError(['Unable to determine number of CT scenarios ' ...
                          'for biological optimization!']);
end

if numOfCtScen < 1
    matRad_cfg.dispError(['Biological optimization requires at least one ' ...
                          'CT scenario!']);
end
end

function values = matRad_asCtScenarioCell(values, numOfCtScen, numOfVoxels, fieldName)
matRad_cfg = MatRad_Config.instance();
if iscell(values)
    values = values(:);
elseif isnumeric(values)
    if numOfCtScen == 1 && isvector(values)
        values = {values(:)};
    elseif size(values, 1) == numOfVoxels && size(values, 2) == numOfCtScen
        values = mat2cell(values, numOfVoxels, ones(1, numOfCtScen));
        values = values(:);
    else
        matRad_cfg.dispError(['%s must be a cell array by CT scenario or ' ...
                              'a numeric array with one column per CT scenario!'], fieldName);
    end
else
    matRad_cfg.dispError(['%s must be numeric or a cell array of numeric ' ...
                          'vectors!'], fieldName);
end
end

function [ax, bx] = matRad_validateAxBxCells(ax, bx, numOfCtScen, numOfVoxels)
matRad_cfg = MatRad_Config.instance();
if ~iscell(ax) || ~iscell(bx) || numel(ax) ~= numOfCtScen || ...
        numel(bx) ~= numOfCtScen
    matRad_cfg.dispError(['dij.ax and dij.bx must have one cell entry per ' ...
                          'CT scenario!']);
end

for scenIx = 1:numOfCtScen
    matRad_validateBiologicalParameterVector(ax{scenIx}, numOfVoxels, ...
                                             sprintf('dij.ax{%d}', scenIx));
    matRad_validateBiologicalParameterVector(bx{scenIx}, numOfVoxels, ...
                                             sprintf('dij.bx{%d}', scenIx));
    ax{scenIx} = ax{scenIx}(:);
    bx{scenIx} = bx{scenIx}(:);
end
end

function matRad_validateBiologicalParameterVector(value, numOfVoxels, fieldName)
matRad_cfg = MatRad_Config.instance();
if ~isnumeric(value) || ~isvector(value) || numel(value) ~= numOfVoxels
    matRad_cfg.dispError(['%s must be a numeric vector with one entry per ' ...
                          'dose-grid voxel!'], fieldName);
end
end

function ixZeroDose = matRad_getZeroDoseVoxelMask(dij, numOfVoxels)
matRad_cfg = MatRad_Config.instance();
if ~isfield(dij, 'physicalDose') || ~iscell(dij.physicalDose)
    matRad_cfg.dispError(['dij.physicalDose is required to prepare ' ...
                          'biological optimization quantities!']);
end

validScen = find(~cellfun(@isempty, dij.physicalDose));
if isempty(validScen)
    matRad_cfg.dispError(['At least one physical dose influence matrix is ' ...
                          'required to prepare biological optimization quantities!']);
end

doseSum = zeros(numOfVoxels, 1);
for scenIx = validScen(:)'
    doseMatrix = dij.physicalDose{scenIx};
    if size(doseMatrix, 1) ~= numOfVoxels
        matRad_cfg.dispError(['dij.physicalDose{%d} has an inconsistent ' ...
                              'number of dose-grid voxels!'], scenIx);
    end
    doseSum = doseSum + doseMatrix * ones(size(doseMatrix, 2), 1);
end
ixZeroDose = doseSum == 0;
end

function matRad_validateAxBxAgainstCst(ax, bx, cst, numOfVoxels)
matRad_cfg = MatRad_Config.instance();
[refAx, refBx] = matRad_getPhotonLQMParameters(cst, numOfVoxels);
if numel(ax) ~= numel(refAx) || numel(bx) ~= numel(refBx)
    matRad_cfg.dispError(['Inconsistent number of biological CT scenarios ' ...
                          'in dij.ax and/or dij.bx!']);
end

isConsistent = cellfun(@(ax1, bx1, ax2, bx2) ...
                       isequal(ax1(ax1 ~= 0), ax2(ax1 ~= 0)) && ...
                       isequal(bx1(bx1 ~= 0), bx2(bx1 ~= 0)), ...
                       ax, bx, refAx, refBx);
if ~all(isConsistent)
    matRad_cfg.dispError(['Inconsistent biological parameters in dij.ax ' ...
                          'and/or dij.bx - please recalculate dose influence matrix before ' ...
                          'optimization!']);
end
end
