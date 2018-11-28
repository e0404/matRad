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

%display options
ompMCoptions.verbose = true;

% start MC control          
ompMCoptions.nHistories = 10000000;
ompMCoptions.nBatches = 10;
ompMCoptions.randomSeeds = [97 33];


%start source definition      
ompMCoptions.spectrumFile = [pwd '/ompMC/spectra/mohan6.spectrum'];
ompMCoptions.monoEnergy = 0.1; 
ompMCoptions.charge = 0;
ompMCoptions.colliBounds = [-2.5 2.5 -2.5 2.5];
ompMCoptions.ssd = 90.0; %This has to be calculated by matRad?
                                                                    
% start MC transport
ompMCoptions.photonXsection = [pwd '/ompMC/data/xcom'];
ompMCoptions.pegsFile = [pwd '/ompMC/pegs4/700icru.pegs4dat'];

ompMCoptions.global_ecut = 0.700;
ompMCoptions.global_pcut = 0.010;        


%% Create Material Density Cube
materialFile = [pwd '/ompMC/700icru.pegs4dat'];
material = cell(4,5);
material{1,1} = 'AIR700ICRU';
material{1,2} = -1024; material{1,3} = -974;
material{1,4} = 0.001; material{1,5} = 0.044;
material{2,1} = 'LUNG700ICRU';
material{2,2} = -974; material{2,3} = -724;
material{2,4} = 0.044; material{2,5} = 0.302;
material{3,1} = 'ICRUTISSUE700ICRU';
material{3,2} = -724; material{3,3} = 101;
material{3,4} = 0.302; material{3,5} = 1.101;
material{4,1} = 'ICRPBONE700ICRU';
material{4,2} = 101; material{4,3} = 1976;
material{4,4} = 1.101; material{4,5} = 2.088;

for s=1:dij.numOfScenarios
        % From HU to densities
    [cubeRho{s},cubeMatIx{s}] = arrayfun(@(HU) HUtoDensityAndMaterial(HU,material),ct.cubeHU{s});
    
    %C Indexing
    cubeMatIx{s} = int32(cubeMatIx{s} - 1);
end

ompMCgeo.material = material;
ompMCgeo.materialFile = materialFile;

scale = 10; %cm?

ompMCgeo.xBounds = [(ct.x - ct.resolution.x*0.5) (ct.x(ct.cubeDim(1)) + ct.resolution.x)] ./ scale;
ompMCgeo.yBounds = [(ct.y - ct.resolution.y*0.5) (ct.y(ct.cubeDim(2)) + ct.resolution.y)] ./ scale;
ompMCgeo.zBounds = [(ct.z - ct.resolution.z*0.5) (ct.z(ct.cubeDim(3)) + ct.resolution.z)] ./ scale;


%% Call the OmpMC interface
%initialize waitbar
figureWait = waitbar(0,'OmpMC photon dose influence matrix calculation..');
fprintf('matRad: OmpMC photon dose calculation... ');

%run over all scenarios
for s = 1:dij.numOfScenarios
    fprintf('Scenario %d... ',s);
     
    ompMCgeo.isoCenter = [stf(:).isoCenter];
    
    dij.physicalDose{s} = matRad_ompInterface(cubeRho{s},cubeMatIx{s},ompMCgeo,ompMCoptions);
end


try
    % wait 0.1s for closing all waitbars
    allWaitBarFigures = findall(0,'type','figure','tag','TMWWaitbar');
    delete(allWaitBarFigures);
    pause(0.1);
catch
end

end

function [rho,matIx] = HUtoDensityAndMaterial(HU,material)
    rho = [];
    for i = 1:size(material,1)
        if HU <= material{i,3}
            matIx = i;
            m = (material{i,5}-material{i,4})/(material{i,3}-material{i,2});
            b = material{i,4} - m*material{i,2};
            rho = m*HU + b;
            break;
        end
    end
    if isempty(rho)
        error(['Error in CT density calculation: No material/density could be assigned to HU value of ' num2str(HU) '!']);
    end
end
