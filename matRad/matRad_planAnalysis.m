function resultGUI = matRad_planAnalysis(resultGUI, ct, cst, stf, pln, varargin)
% matRad_planAnalysis calculates DVH and quality indicators for a plan.
%
% call
%   resultGUI = matRad_planAnalysis(resultGUI,ct,cst,stf,pln)
%   resultGUI = matRad_planAnalysis(resultGUI,ct,cst,stf,pln,varargin)
%
% input
%   resultGUI:              matRad resultGUI struct containing dose cubes
%   ct:                     matRad ct struct with computed tomography data
%   cst:                    matRad cst cell array with structure definitions
%   stf:                    matRad stf struct with beam information
%   pln:                    matRad pln struct with plan information
%   name / value pairs:     Optional parameters for analysis customization
%   quantity: (optional)    resultGUI dose quantity to analyse
%   evaluationMode:(optional) 'perFraction' or 'total' for figures/tables
%   doseWindow: (optional)  dose axis window for DVH display
%   refGy: (optional)       Per-fraction dose values for V_XGy calculation
%   refVol:(optional)       Volume percentages for D_X calculation
%
% output
%   resultGUI:              Updated resultGUI with analysis data
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2024-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser();
p.addRequired('ct', @isstruct);
p.addRequired('cst', @iscell);
p.addRequired('stf', @isstruct);
p.addRequired('pln', @isstruct);
p.addParameter('refGy', [], @isnumeric);
p.addParameter('refVol', [2 5 95 98], @isnumeric);
p.addParameter('showDVH', true, @islogical);
p.addParameter('showQI', true, @islogical);
p.addParameter('quantity', '', @(x) ischar(x) || isstring(x));
p.addParameter('evaluationMode', 'perFraction', @(x) ischar(x) || isstring(x));
p.addParameter('doseWindow', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
p.parse(ct, cst, stf, pln, varargin{:});

ct = p.Results.ct;
cst = p.Results.cst;
pln = p.Results.pln;
refGy = p.Results.refGy;
refVol = p.Results.refVol;
showDVH = p.Results.showDVH;
showQI = p.Results.showQI;
quantity = p.Results.quantity;
doseWindow = p.Results.doseWindow;
evaluationMode = p.Results.evaluationMode;

[~, evaluationMode, evaluationScale] = matRad_convertToEvaluationMode([], pln, evaluationMode);

if ~isempty(doseWindow)
    doseWindow = doseWindow(:)';
end

visQ = matRad_resolveDoseAnalysisQuantity(resultGUI, pln, quantity);

if ~isfield(pln, 'multScen')
    pln.multScen = 'nomScen';
end

if ~isa(pln.multScen, 'matRad_ScenarioModel')
    pln.multScen = matRad_ScenarioModel.create(pln.multScen, ct);
end

doseCube = resultGUI.(visQ);
resultGUI.dvh = matRad_calcDVH(cst, doseCube, 'cum');
resultGUI.qi  = matRad_calcQualityIndicators(cst, pln, doseCube, refGy, refVol);
resultGUI.analysisQuantity = visQ;
resultGUI.evaluationModeBase = 'perFraction';
resultGUI.evaluationMode = evaluationMode;
resultGUI.evaluationScale = evaluationScale;
resultGUI.displayDvh = matRad_convertDvhForEvaluation(resultGUI.dvh, pln, evaluationMode);
resultGUI.displayQi = matRad_convertQiForEvaluation(resultGUI.qi, pln, evaluationMode);

dvhScen = {};
numScenarios = pln.multScen.totNumScen;
if numScenarios > 1
    for i = 1:numScenarios
        scenFieldName = sprintf('%s_scen%d', visQ, i);
        if isfield(resultGUI, scenFieldName)
            dvhScen{i} = matRad_convertDvhForEvaluation( ...
                                                        matRad_calcDVH(cst, resultGUI.(scenFieldName), 'cum'), pln, evaluationMode); %#ok<AGROW>
        end
    end
end

if showDVH || showQI
    matRad_cfg = MatRad_Config.instance();
    figure('Color', matRad_cfg.gui.backgroundColor);

    colorSpec = {'Color', matRad_cfg.gui.elementColor, ...
                 'XColor', matRad_cfg.gui.textColor, ...
                 'YColor', matRad_cfg.gui.textColor, ...
                 'GridColor', matRad_cfg.gui.textColor, ...
                 'MinorGridColor', matRad_cfg.gui.backgroundColor};

    if showDVH && showQI
        hDVHax = subplot(2, 1, 1, colorSpec{:});
        hQIax = subplot(2, 1, 2, colorSpec{:});
    elseif showDVH
        hDVHax = subplot(1, 1, 1, colorSpec{:});
    elseif showQI
        hQIax = subplot(1, 1, 1, colorSpec{:});
    end
end

if showDVH
    matRad_showDVH(resultGUI.displayDvh, cst, pln, 'axesHandle', hDVHax, 'LineWidth', 3);

    for i = 1:numel(dvhScen)
        matRad_showDVH(dvhScen{i}, cst, pln, 'axesHandle', hDVHax, ...
                       'LineWidth', 0.5, 'plotLegend', false, 'LineStyle', '--');
    end

    if ~isempty(doseWindow)
        xlim(hDVHax, doseWindow);
    end
end

if showQI
    matRad_showQualityIndicators(hQIax, resultGUI.displayQi);
end

end

function dvh = matRad_convertDvhForEvaluation(dvh, pln, evaluationMode)
if isempty(dvh)
    return
end

for i = 1:numel(dvh)
    if isfield(dvh(i), 'doseGrid') && ~isempty(dvh(i).doseGrid)
        dvh(i).doseGrid = matRad_convertToEvaluationMode( ...
                                                         dvh(i).doseGrid, pln, evaluationMode);
    end
end
end

function qi = matRad_convertQiForEvaluation(qi, pln, evaluationMode)
if isempty(qi)
    return
end

doseFields = {'mean', 'std', 'max', 'min', 'referenceDose'};
for i = 1:numel(qi)
    for j = 1:numel(doseFields)
        fieldName = doseFields{j};
        if isfield(qi(i), fieldName) && isnumeric(qi(i).(fieldName))
            qi(i).(fieldName) = matRad_convertToEvaluationMode( ...
                                                               qi(i).(fieldName), pln, evaluationMode);
        end
    end

    fields = fieldnames(qi(i));
    for j = 1:numel(fields)
        fieldName = fields{j};
        if strncmp(fieldName, 'D_', 2) && isnumeric(qi(i).(fieldName))
            qi(i).(fieldName) = matRad_convertToEvaluationMode( ...
                                                               qi(i).(fieldName), pln, evaluationMode);
        end
    end
end
end
