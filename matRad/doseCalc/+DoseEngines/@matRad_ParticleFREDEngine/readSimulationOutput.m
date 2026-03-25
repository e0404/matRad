function [doseCube, letCube, loadFileName] = readSimulationOutput(runFolder,calcDoseDirect,calcLET)
% FRED helper to read simulation output
% call:
%   readSimulationOutput(runFolder,calcDoseDirect, varargin)
% 
% input:
%   runFolder:          path to folder containing the simulation files
%   calcDoseDirect:     boolean to trigger dij or .mhd reading
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matRad_cfg = MatRad_Config.instance();

p = inputParser();
addRequired(p, 'runFolder', @ischar);
addRequired(p, 'calcDoseDirect', @islogical);
addOptional(p, 'calcLET',false,@islogical);

parse(p, runFolder,calcDoseDirect, calcLET);

runFolder       = p.Results.runFolder;
calcDoseDirect  = p.Results.calcDoseDirect;
calcLET         = p.Results.calcLET;

doseCube = [];
letCube  = [];

outFolder = fullfile(runFolder,'out');

if ~calcDoseDirect
    letCompatible = true;

    scoreIjFolder = fullfile(outFolder, 'scoreij');
    doseDijFile = 'Phantom.Dose.bin';
    letdDijFile = 'Phantom.LETd.bin';
    if ~exist(scoreIjFolder,"dir")
        scoreIjFolder = outFolder;
        doseDijFile = 'Dij.bin';
        letCompatible = false;
    end
    loadFileName = fullfile(scoreIjFolder,doseDijFile);
    
    matRad_cfg.dispInfo('Looking for scorer-ij output in sub folder: %s\n', strrep(scoreIjFolder, '\', '\\'));
    
    % read dij matrix
    if isfile(loadFileName)
        dijMatrices = DoseEngines.matRad_ParticleFREDEngine.readSparseDijBin(loadFileName);
        doseCube = dijMatrices{1};
    else
        matRad_cfg.dispError(sprintf('Unable to find file: %s', strrep(loadFileName, '\', '\\')));
    end

    if calcLET && letCompatible

        
        letdDijFileName = fullfile(scoreIjFolder,letdDijFile);
        
        try
            dijMatrices = DoseEngines.matRad_ParticleFREDEngine.readSparseDijBin(letdDijFileName);
            % TODO: store nominator and denominator in dij and allow this
            % structure for LETd calculation in calcCubes and
            % backProjection.
            dijNom = dijMatrices{1};
            dijDen = dijMatrices{2};
            letCube = spfun(@(x) 1./x,dijDen);
            letCube = letCube .* dijNom ./ 10; %divided by 10 to have keV/um
        catch
            matRad_cfg.dispError('unable to load file: %s',letdDijFileName);
        end
    end
else
    scoreFolder = fullfile(outFolder, 'score');
    doseCubeFileName = 'Phantom.Dose.mhd';
    letdCubeFileName = 'Phantom.LETd.mhd';
    if ~exist(scoreFolder,"dir")
        scoreFolder = outFolder;
        doseCubeFileName = 'Dose.mhd';
        letdCubeFileName = 'LETd.mhd';
    end
    
    loadFileName = fullfile(scoreFolder, doseCubeFileName);

    matRad_cfg.dispInfo(sprintf('Looking for scorer file in sub folder: %s\n', strrep(scoreFolder, '\', '\\')));

    if isfile(loadFileName)
        doseCube = matRad_readMHD(loadFileName);
    else
        matRad_cfg.dispError(sprintf('Unable to find file: %s', strrep(loadFileName, '\', '\\')));
    end

    if calcLET
        try 
            letCube = matRad_readMHD(fullfile(scoreFolder, letdCubeFileName));
        catch
            matRad_cfg.dispError('Unable to load file: %s',fullfile(letdDijFolder, letdCubeFileName));
        end
    end
    
end

matRad_cfg.dispInfo('Loading succesful!\n');
end
