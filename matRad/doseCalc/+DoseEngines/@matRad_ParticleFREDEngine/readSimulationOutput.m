function [doseCube, letCube] = readSimulationOutput(runFolder,calcDoseDirect, varargin)

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

if ~calcDoseDirect

    doseDijFolder = fullfile(runFolder, 'out', 'scoreij');
    doseDijFile = 'Phantom.Dose.bin';
    
    doseDijFileName = fullfile(doseDijFolder,doseDijFile);
    
    % read dij matrix
 
    doseCube = DoseEngines.matRad_ParticleFREDEngine.readSparseDijBin(doseDijFileName);

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
    
    doseCube = matRad_readMHD(fullfile(doseCubeFolder, doseCubeFileName));
        
    if calcLET

        letdDijFolder   = doseCubeFolder;
        letdCubeFileName     = 'Phantom.LETd.mhd';

        try 
            letCube = matRad_readMHD(fullfile(letdDijFolder, letdCubeFileName));
        catch
            matRad_cfg.dispError('unable to load file: %s',fullfile(letdDijFolder, letdCubeFileName));
        end
    end
    
end
end