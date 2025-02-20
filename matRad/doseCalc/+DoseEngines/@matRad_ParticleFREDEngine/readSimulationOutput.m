function [doseCube, letCube, loadFileName] = readSimulationOutput(runFolder,calcDoseDirect,varargin)
% FRED helper to read simulation output
% call
%   readSimulationOutput(runFolder,calcDoseDirect, varargin)
% 
% input
%   runFolder:          path to folder containing the simulation files
%   calcDoseDirect:     boolean to trigger dij or .mhd reading
%   
%  optional:
%   calLET:             addirional boolean to trigger loading of LETd
%   readFunctionHandle: handle to readout function for ij-scorer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team.
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
addParameter(p, 'calcLET',0,@islogical);
addParameter(p, 'readFunctionHandle', @(lFile) DoseEngines.matRad_ParticleFREDEngine.readSparseDijBin(lFile))

parse(p, runFolder,calcDoseDirect, varargin{:});

runFolder       = p.Results.runFolder;
calcDoseDirect  = p.Results.calcDoseDirect;
calcLET         = p.Results.calcLET;

doseCube = [];
letCube  = [];
loadFileName = [];

if ~calcDoseDirect

    doseDijFolder = fullfile(runFolder, 'out', 'scoreij');
    doseDijFile = 'Phantom.Dose.bin';
    loadFileName = fullfile(doseDijFolder,doseDijFile);
    
    matRad_cfg.dispInfo(sprintf('Looking for scorer-ij output in sub folder: %s\n', strrep(doseDijFolder, '\', '\\')));
    
    % read dij matrix
    if isfile(loadFileName)
%        doseCube = DoseEngines.matRad_ParticleFREDEngine.readSparseDijBin(loadFileName);
        doseCube = p.Results.readFunctionHandle(loadFileName);
    else
        matRad_cfg.dispError(sprintf('Unable to find file: %s', strrep(loadFileName, '\', '\\')));
    end

    if calcLET
        
        letdDijFile = 'Phantom.LETd.bin';
        letdDijFileName = fullfile(doseDijFolder,letdDijFile);
        
        try
%            letCube = DoseEngines.matRad_ParticleFREDEngine.readSparseDijBin(letdDijFileName);
            letCube = p.Results.readFunctionHandle(letdDijFileName);
        catch
            matRad_cfg.dispError('unable to load file: %s',letdDijFileName);
        end
    end
else
    
    doseCubeFolder = fullfile(runFolder, 'out', 'score');
    doseCubeFileName = 'Phantom.Dose.mhd';
    loadFileName = fullfile(doseCubeFolder, doseCubeFileName);

    matRad_cfg.dispInfo(sprintf('Looking for scorer file in sub folder: %s\n', strrep(doseCubeFolder, '\', '\\')));

    if isfile(loadFileName)
        doseCube = matRad_readMHD(loadFileName);
    else
        matRad_cfg.dispError(sprintf('Unable to find file: %s', strrep(loadFileName, '\', '\\')));
    end

    if calcLET

        letdDijFolder   = doseCubeFolder;
        letdCubeFileName     = 'Phantom.LETd.mhd';

        try 
            letCube = matRad_readMHD(fullfile(letdDijFolder, letdCubeFileName));
        catch
            matRad_cfg.dispError('Unable to load file: %s',fullfile(letdDijFolder, letdCubeFileName));
        end
    end
    
end

matRad_cfg.dispInfo('Loading succesful!\n');
end