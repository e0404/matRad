function resultGUI = matRad_planAnalysis(resultGUI,ct,cst,stf,pln,varargin)
% matRad plan analysis function
% This function performs analysis on radiation therapy plans, including DVH (Dose-Volume Histogram) and quality indicators.
% It optionally displays these analyses based on input parameters.
%
% input
%   resultGUI:              matRad resultGUI struct containing the analysis results
%   ct:                     matRad ct struct with computed tomography data
%   cst:                    matRad cst cell array with structure definitions
%   stf:                    matRad stf struct with beam information
%   pln:                    matRad pln struct with plan information
%   name / value pairs:     Optional parameters for analysis customization
%   refGy: (optional)       Dose values for V_XGy calculation (default: [40 50 60])
%   refVol:(optional)       Volume percentages for D_X calculation (default: [2 5 95 98])
%
% output
%   resultGUI:              Updated resultGUI with analysis data

% Initialize input parser for function arguments
p = inputParser();

% Define required inputs
p.addRequired('ct',@isstruct);
p.addRequired('cst',@iscell);
p.addRequired('stf',@isstruct);
p.addRequired('pln',@isstruct);

%
% Copyright 2024 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialize input parser for optional parameters
p = inputParser();

% Define required inputs again for clarity
p.addRequired('ct',@isstruct);
p.addRequired('cst',@iscell);
p.addRequired('stf',@isstruct);
p.addRequired('pln',@isstruct);

% Define optional parameters with default values
p.addParameter('refGy',[40 50 60],@isnumeric); % Reference dose values for V_XGy calculation
p.addParameter('refVol',[2 5 95 98],@isnumeric); % Reference volume percentages for D_X calculation
p.addParameter('showDVH',true,@islogical); % Flag to show or hide the DVH plot
p.addParameter('showQI',true,@islogical); % Flag to show or hide the Quality Indicators plot

% Parse input arguments to extract values
p.parse(ct,cst,stf,pln,varargin{:});

% Assign parsed values to variables
ct = p.Results.ct;
cst = p.Results.cst;
stf = p.Results.stf;
pln = p.Results.pln;
refGy = p.Results.refGy;
refVol = p.Results.refVol;
showDVH = p.Results.showDVH;
showQI = p.Results.showQI;

% Determine which dose cube to use based on resultGUI structure
if isfield(resultGUI,'RBExD')
    doseCube = resultGUI.RBExD; % Use RBExD if available
else
    doseCube = resultGUI.physicalDose; % Default to physicalDose otherwise
end

% Calculate DVH and quality indicators
resultGUI.dvh = matRad_calcDVH(cst,doseCube,'cum'); % Calculate cumulative DVH
resultGUI.qi  = matRad_calcQualityIndicators(cst,pln,doseCube,refGy,refVol); % Calculate quality indicators

% Configuration for GUI appearance
matRad_cfg = MatRad_Config.instance();

% Create figure for plots with background color from configuration
hF = figure('Color',matRad_cfg.gui.backgroundColor);

% Determine subplot layout based on flags
if showDVH && showQI
    hDVHax = subplot(2,1,1); % DVH plot area
    hQIax = subplot(2,1,2); % Quality Indicators plot area
elseif showDVH
    hDVHax = subplot(1,1,1); % Only DVH plot
elseif showQI
    hQIax = subplot(1,1,1); % Only Quality Indicators plot
end

% Display DVH if enabled
if showDVH
    matRad_showDVH(hDVHax,resultGUI.dvh,cst,pln); % Show DVH plot
end

% Display Quality Indicators if enabled
if showQI
    matRad_showQualityIndicators(hQIax,resultGUI.qi); % Show Quality Indicators plot
end

end