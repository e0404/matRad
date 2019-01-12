function dij = matRad_calcParticleDoseMC2(ct,stf,pln,cst)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad particle dose calculation wrapper using MCsquare
%
% call
%   dij = matRad_calcParticleDoseMC2(ct,stf,pln,cst)
%
% input
%   ct:             ct cube
%   stf:            matRad steering information struct
%   pln:            matRad plan meta information struct
%   cst:            matRad cst struct
%
% output
%   dij:            matRad dij struct
%
% References
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("Dose influence matrix computation using MCsquare")
% Copy MCsquare code locally
copyfile([pwd '/MC2/*'],'.');

cleanUp = false;

nbThreads = 4;
nbParallel = 8;%nproc();
if nbParallel > 4
  nbParallel = 2*ceil(nbParallel/nbThreads);
end

% meta information for dij
dij.numOfBeams         = pln.propStf.numOfBeams;
dij.numOfVoxels        = prod(ct.cubeDim);
dij.resolution         = ct.resolution;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.dimensions         = ct.cubeDim;
dij.numOfScenarios     = 1;
dij.optimExtra.size    = 0;
dij.optimExtra.quantity = {};

% set consistent random seed (enables reproducibility)
%rng(0);
rand('state',0)

% set number of particles simulated per pencil beam
nCasePerBixel              = 1000;

% set relative dose cutoff for storage in dose influence matrix
relDoseCutoff = 10^(-3);
% set absolute calibration factor
% convert from sum of dose in Gy for all histories to Gy per 1e6 primaries
absCalibrationFactorMC2 = 1.0e+6/nCasePerBixel;

% set up arrays for book keeping
dij.bixelNum = NaN*ones(dij.totalNumOfRays,1);
dij.rayNum   = NaN*ones(dij.totalNumOfRays,1);
dij.beamNum  = NaN*ones(dij.totalNumOfRays,1);

% Allocate space for dij.physicalDose sparse matrix
for i = 1:dij.numOfScenarios
    dij.physicalDose{i} = spalloc(prod(ct.cubeDim),dij.totalNumOfBixels,1);
end

% helper function for energy selection
round2 = @(a,b)round(a*10^b)/10^b;

% Allocate memory for dose_temp cell array
numOfBixelsContainer = nbParallel;

if ~strcmp(pln.radiationMode,'protons')
    errordlg('MCsquare is only supported for protons');
end

doseTmpContainer = cell(numOfBixelsContainer,dij.numOfScenarios);
if isequal(pln.propOpt.bioOptimization,'const_RBExD')
            dij.RBE = 1.1;
            fprintf(['matRad: Using a constant RBE of 1.1 \n']);
end

% Only take voxels inside patient.
V = [cst{:,4}];
V = unique(vertcat(V{:}));

% load machine file
fileName = [pln.radiationMode '_' pln.machine];
try
   load([fileparts(mfilename('fullpath')) filesep fileName]);
catch
   error(['Could not find the following machine file: ' fileName ]);
end

%%%%%%%%%%%%%%%%%%%%%%

% compute SSDs
stf = matRad_computeSSD(stf,ct);

fprintf('matRad: Simply book Dij and creates MCsquare input files...\n');
counter = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MCbixels = [];
MCparameters.nbFields = dij.numOfBeams;
MCparameters.nbThreads = nbThreads;

for i = 1:dij.numOfBeams % loop over all beams

    fprintf(['Beam ' num2str(i) ' of ' num2str(dij.numOfBeams) ': \n']);

    bixelsPerBeam = 0;

    for j = 1:stf(i).numOfRays % loop over all rays

        if ~isempty(stf(i).ray(j).energy) % ray energy check

            for k = 1:stf(i).numOfBixelsPerRay(j) % loop over all bixels per ray

                counter = counter + 1;
                bixelsPerBeam = bixelsPerBeam + 1;

                % remember beam and bixel number
                dij.beamNum(counter)  = i;
                dij.rayNum(counter)   = j;
                dij.bixelNum(counter) = k;

                % find energy index in base data
                energyIx = find(round2(stf(i).ray(j).energy(k),4) == round2([machine.data.energy],4));
                focusIx=stf(i).ray(j).focusIx(k);
                rayPos_bev = stf(i).ray(j).rayPos_bev;

                posX = -rayPos_bev(3);
                posY = rayPos_bev(1);

                % MC control
                MCbixel = [];

                MCbixel.ncase  = nCasePerBixel; % number of histories

                % create a random seed for every bixel
                MCbixel.rngSeed = randi(50000);
                MCbixel.pln.beamNb  = i;
                MCbixel.pln.rayNb   = j;
                MCbixel.pln.bixelNb = k;
                MCbixel.pln.counter = counter;
                MCbixel.stf.posX   = posX;
                MCbixel.stf.posY   = posY;
                MCbixel.stf.energy = stf(i).ray(j).energy(k);
                MCbixel.stf.focus  = machine.data(energyIx).initFocus.SisFWHMAtIso(focusIx);

                MCbixels = [MCbixels,MCbixel];
            end % loop over all bixels per ray
        end % ray energy check
    end % loop over all rays
end % loop over all beams

if counter ~= dij.totalNumOfBixels
  error('counter != dij.totalNumOfBixels');
end

MCparameters.bixels = MCbixels; 
% save data for interface to MCsquare
%if exist('OCTAVE_VERSION','builtin');
  % OCTAVE
  save('-v7','matRad_data.mat','ct','pln','MCparameters');

%else
  % MATLAB
  %save('matRad_data.mat','ct','pln','MCparameters');
%end
system(['C:\Users\localbangertm\Anaconda3\python.exe ' pwd '/calc_pencilBeamMCsquare.py --exportCT']);

fprintf('matRad: Particle dose calculation using MCsquare...\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for counter = 1:numOfBixelsContainer:dij.totalNumOfBixels

    nbRuns = numOfBixelsContainer;
    if (counter+nbRuns)> dij.totalNumOfBixels
      nbRuns = dij.totalNumOfBixels-counter+1;
    end

    fprintf(['Running chunk from ' num2str(counter) ' to ' num2str(counter+nbRuns-1) ' out of ' num2str(dij.totalNumOfBixels) '\n']);

    if exist('OCTAVE_VERSION') ~= 0
        % pkg install -forge parallel
        pkg load parallel
        inputs = counter:counter+nbRuns-1;
        resultsArray = [];
        [resultsArray] = pararrayfun (nbParallel, @(runNb) computePencilDose(runNb,cleanUp,ct.cubeDim,dij.numOfVoxels,V,relDoseCutoff,absCalibrationFactorMC2), inputs);
         % Save dose for every bixel in cell array
        doseTmpContainer = reshape(mat2cell(resultsArray,[rows(resultsArray)],ones(1,nbRuns)),[nbRuns 1]);
    else
      %parfor runNb = counter:counter+nbRuns-1
      for runNb = counter:counter+nbRuns-1
          % Save dose for every bixel in cell array
          doseTmpContainer{mod(runNb-1,numOfBixelsContainer)+1,1} = computePencilDose(runNb,cleanUp,ct.cubeDim,dij.numOfVoxels,V,relDoseCutoff,absCalibrationFactorMC2);
      end
    end
     
    % save computation time and memory by sequentially filling the
    % sparse matrix dose.dij from the cell array
    dij.physicalDose{1}(:,counter:counter+nbRuns-1) = [doseTmpContainer{1:nbRuns,1}];

    % Display progress and update text only 200 times
    if mod(counter+nbRuns-1,max(1,round(dij.totalNumOfBixels/200))) == 0
        matRad_progress((counter+nbRuns-1)/max(1,round(dij.totalNumOfBixels/200)),...
                        floor(dij.totalNumOfBixels/max(1,round(dij.totalNumOfBixels/200))));
        fprintf('\n');          
    end

end

end

function sparseDose = computePencilDose(runNb,cleanUp,cubeDim,numOfVoxels,V,relDoseCutoff,absCalibrationFactorMC2)

    % MCsquare beamlet input data
    [status, cmdout] = system(['C:\Users\localbangertm\Anaconda3\python.exe ' pwd '/calc_pencilBeamMCsquare.py --bixelNb ' num2str(runNb)]);
    if status ~= 0
        error(['Simulation did not return success:' cmdout]);
    end

    data = load(['bixelDose_' num2str(runNb) '.mat']);
    bixelDose = data.bixelDose;

    % MCsquare import dose cube
    % Clean-up data
    if cleanUp
        %system(['rm MCpencilbeam_bixelNb_' num2str(runNb) '.txt']);
        %system(['rm config_bixelNb_' num2str(runNb) '.txt']);
        %system(['rm -r Outputs_bixelNb_' num2str(runNb)]);
        %system(['rm bixelDose_' num2str(runNb) '.mat']);
        rmdir(['Outputs_bixelNb_' num2str(runNb)],'s');
        delete(['MCpencilbeam_bixelNb_' num2str(runNb) '.txt']);
        delete(['config_bixelNb_' num2str(runNb) '.txt']);
        delete(['bixelDose_' num2str(runNb) '.mat']);
    end

    % apply relative dose cutoff
    doseCutoff                        = relDoseCutoff*max(bixelDose(:));
    bixelDose(bixelDose < doseCutoff) = 0;

    % apply absolute calibration factor
    bixelDose = bixelDose*absCalibrationFactorMC2;

    sparseDose = sparse(V,1,double(bixelDose(V)),numOfVoxels,1);
end

