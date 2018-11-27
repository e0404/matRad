function dij = matRad_calcPhotonDoseOmpMC(ct,stf,pln,cst)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad ompMC monte carlo photon dose calculation wrapper
% 
% call
%   dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst)
%
% input
%   ct:                         matRad ct struct
%   stf:                        matRad steering information struct
%   pln:                        matRad plan meta information struct
%   cst:                        matRad cst struct
% output
%   dij:                        matRad dij struct
%
% References
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default: dose influence matrix computation
if ~exist('calcDoseDirect','var')
    calcDoseDirect = false;
end

% set output level. 0 = no vmc specific output. 1 = print to matlab cmd.
% 2 = open in terminal(s)
verbose = 1;

if ~isdeployed % only if _not_ running as standalone    
    % add path for optimization functions
    matRadRootDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(matRadRootDir,'ompMC'))
end

% meta information for dij
dij.numOfBeams         = pln.propStf.numOfBeams;
dij.numOfVoxels        = prod(ct.cubeDim);
dij.resolution         = ct.resolution;
dij.dimensions         = ct.cubeDim;
dij.numOfScenarios     = 1;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);

% check if full dose influence data is required
%{
if calcDoseDirect 
    numOfColumnsDij           = length(stf);
    numOfBixelsContainer = 1;
else
    numOfColumnsDij           = dij.totalNumOfBixels;
    numOfBixelsContainer = ceil(dij.totalNumOfBixels/10);
end
%}

numOfColumnsDij           = dij.totalNumOfBixels;
numOfBixelsContainer = ceil(dij.totalNumOfBixels/10);

% set up arrays for book keeping
dij.bixelNum = NaN*ones(numOfColumnsDij,1);
dij.rayNum   = NaN*ones(numOfColumnsDij,1);
dij.beamNum  = NaN*ones(numOfColumnsDij,1);

bixelNum = NaN*ones(dij.totalNumOfBixels,1);
rayNum   = NaN*ones(dij.totalNumOfBixels,1);
beamNum  = NaN*ones(dij.totalNumOfBixels,1);

doseTmpContainer = cell(numOfBixelsContainer,dij.numOfScenarios);

% Allocate space for dij.physicalDose sparse matrix
%for i = 1:dij.numOfScenarios
%    dij.physicalDose{i} = spalloc(prod(ct.cubeDim),numOfColumnsDij,1);
%end



% set environment variables for ompMC
%{
if exist(['vmc++' filesep 'bin'],'dir') ~= 7
    error(['Could not locate vmc++ environment. ' ...
          'Please provide the files in the correct folder structure at matRadroot' filesep 'vmc++.']);
else
    VMCPath     = fullfile(pwd , 'vmc++');
    runsPath    = fullfile(VMCPath, 'runs');
    phantomPath = fullfile(VMCPath, 'phantoms');

    setenv('vmc_home',VMCPath);
    setenv('vmc_dir',runsPath);
    setenv('xvmc_dir',VMCPath);
    
    if isunix
        system(['chmod a+x ' VMCPath filesep 'bin' filesep 'vmc_Linux.exe']);
    end
    
end
%}

%% Setup OmpMC options / parameters
ompMCgeo.nFields = numel(stf);
ompMCgeo.nBixels = [stf(:).totalNumOfBixels];
ompMCgeo.ctRes = [ct.resolution.x ct.resolution.y ct.resolution.z];


ompMCoptions.nHistories = 5000;
ompMCoptions.nBatches = 10;
%ompMCoptions.randomSeed = 0;

ompMCoptions.verbose = true;

%% Call the OmpMC interface
%initialize waitbar
figureWait = waitbar(0,'OmpMC photon dose influence matrix calculation..');
fprintf('matRad: OmpMC photon dose calculation... ');

%run over all scenarios
for s = 1:dij.numOfScenarios
    fprintf('Scenario %d... ',s);
    
    ompMCgeo.isoCenter = [stf(:).isoCenter];
    
    dij.physicalDose{s} = matRad_ompInterface(ct.cube{s},ompMCgeo,ompMCoptions);
end


try
  % wait 0.1s for closing all waitbars
  allWaitBarFigures = findall(0,'type','figure','tag','TMWWaitbar'); 
  delete(allWaitBarFigures);
  pause(0.1); 
catch
end
