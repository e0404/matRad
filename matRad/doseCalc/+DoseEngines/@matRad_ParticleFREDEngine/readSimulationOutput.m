function [doseCube, letCube, loadFileName] = readSimulationOutput(runFolder,calcDoseDirect, varargin)

matRad_cfg = MatRad_Config.instance();

p = inputParser();
addRequired(p, 'runFolder', @ischar);
addRequired(p, 'calcDoseDirect', @islogical);
addOptional(p, 'calcLET',0,@islogical);

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
        doseCube = DoseEngines.matRad_ParticleFREDEngine.readSparseDijBin(loadFileName);
    else
        matRad_cfg.dispError(sprintf('Unable to find file: %s', strrep(loadFileName, '\', '\\')));
    end

    if calcLET
        
        letdDijFile = 'Phantom.LETd.bin';
        letdDijFileName = fullfile(doseDijFolder,letdDijFile);
        
        try
            letCube = DoseEngines.matRad_ParticleFREDEngine.readSparseDijBin(letdDijFileName);
        catch
            matRad_cfg.dispError('unable to load file: %s',letdDijFileName);
        end
    end
else
    
    doseCubeFolder = fullfile(runFolder, 'out', 'score');
    doseCubeFileName = 'Phantom.Dose.mhd';
    loadFileName = fullfile(doseCubeFolder, doseCubeFileName);

    matRad_cfg.dispInfo(sprintf('Looking for scorer-ij file in sub folder: %s\n', strrep(doseCubeFolder, '\', '\\')));

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