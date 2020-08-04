function dij = matRad_calcParticleDoseMCtopas(ct,stf,pln,cst,nCasePerBixel,calcDoseDirect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad TOPAS Monte Carlo proton dose calculation wrapper
%
% call
%   dij = matRad_calcParticleDoseMCtopas(ct,stf,pln,cst,nCasePerBixel,calcDoseDirect)
%
% input
%   ct:                         matRad ct struct
%   stf:                        matRad steering information struct
%   pln:                        matRad plan meta information struct
%   cst:                        matRad cst struct
%   nCasePerBixel               number of histories per beamlet (nCasePerBixel > 1),
%                               max stat uncertainity (0 < nCasePerBixel < 1)
%   calcDoseDirect:             binary switch to enable forward dose
%                               calcualtion
% output
%   dij:                        matRad dij struct
%
% References
%
%   https://aapm.onlinelibrary.wiley.com/doi/abs/10.1118/1.4943377
%   http://www.openmcsquare.org/
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global matRad_cfg;
matRad_cfg = MatRad_Config.instance();

if nargin < 5
    % set number of particles simulated per pencil beam
    nCasePerBixel = matRad_cfg.propMC.MCsquare_defaultHistories;
    matRad_cfg.dispInfo('Using default number of Histories per Bixel: %d\n',nCasePerBixel);
end

if nargin < 6
    calcDoseDirect = false;
end

if isfield(pln,'propMC') && isfield(pln.propMC,'outputVariance')
    matRad_cfg.dispWarning('Variance scoring for TOPAS not yet supported.');
end

if ~calcDoseDirect
    matRad_cfg.dispError('matRad so far only supports direct dose calculation for TOPAS!\n');
end

env = matRad_getEnvironment();

%% Initialize dose Grid as usual
matRad_calcDoseInit;

% for TOPAS we explicitly downsample the ct to the dose grid (might not
% be necessary in future versions with separated grids)
for s = 1:dij.numOfScenarios
    ctResampled = ct;
    ctResampled.cubeHU{s} =  matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y',  dij.ctGrid.z,ct.cubeHU{s}, ...
                                dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'linear');
    ctResampled.resolution = dij.doseGrid.resolution;
    ctResampled.cubeDim = dij.doseGrid.dimensions;
    ctResampled.x = dij.doseGrid.x;
    ctResampled.y = dij.doseGrid.y;
    ctResampled.z = dij.doseGrid.z;
end

%% sending data to topas

load([pln.radiationMode,'_',pln.machine]);
topasConfig = MatRad_TopasConfig();
topasBaseData = MatRad_TopasBaseData(machine,stf);%,TopasConfig);

%Collect weights
w = zeros(sum([stf(:).totalNumOfBixels]),1);
ct = 1;
for i = 1:length(stf)
    for j = 1:stf(i).numOfRays
        rayBix = stf(i).numOfBixelsPerRay(j);
        w(ct:ct+rayBix-1) = stf(i).ray(j).weight;
        ct = ct + rayBix;
    end
end

topasConfig.numHistories = nCasePerBixel;
topasConfig.numOfRuns = matRad_cfg.propMC.topas_defaultNumBatches;
topasConfig.writeAllFiles(ctResampled,pln,stf,topasBaseData,w);

%topasConfig.parallelRuns = true;
%topasConfig.numThreads = 20 / topasConfig.numOfRuns;

%Run Simulation
currDir = cd;
cd(topasConfig.workingDir);
for beamIx = 1:numel(stf)
    for runIx = 1:topasConfig.numOfRuns       
        fname = sprintf('%s_field%d_run%d',topasConfig.label,beamIx,runIx);
        topasCall = sprintf('%s %s.txt > %s.out 2> %s.err',topasConfig.topasExecCommand,fname,fname,fname);
        if topasConfig.parallelRuns
            finishedFiles{runIx} = sprintf('%s.finished',fname);
            delete(finishedFiles{runIx});
            topasCall = [topasCall '; touch ' finishedFiles{runIx} ' &'];
        end
        matRad_cfg.dispInfo('Calling TOPAS: %s\n',topasCall);
        [status,cmdout] = system(topasCall,'-echo');
        if status == 0
            matRad_cfg.dispInfo('TOPAS simulation completed succesfully\n');
        else
            matRad_cfg.dispError('TOPAS simulation exited with error code %d\n',status);
        end
    end
    
    if topasConfig.parallelRuns
        runsFinished = false;
        pause('on');
        while ~runsFinished
            pause(1);
            runsFinished = true;
            fin = cellfun(@(f) exist(f,'file'),finishedFiles);
            runsFinished = all(fin);
        end
    end
        
end
cd(currDir);

%% read out volume scorers from topas simulation
topasCubes = matRad_readTopasData(topasConfig.workingDir);

fnames = fieldnames(topasCubes);
dij.MC_tallies = fnames;
for f = 1:numel(fnames)
    dij.(fnames{f}){1} = topasCubes.(fnames{f});    
end
