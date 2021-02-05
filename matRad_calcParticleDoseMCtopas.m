function dij = matRad_calcParticleDoseMCtopas(ct,stf,pln,cst,nCasePerBixel,calcDoseDirect,openStack)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad TOPAS Monte Carlo proton dose calculation wrapper
%   This calls a TOPAS installation (not included in matRad due to
%   licensing model of TOPAS) for MC simulation
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

matRad_cfg = MatRad_Config.instance();

if nargin < 5
    % set number of particles simulated per pencil beam
    nCasePerBixel = matRad_cfg.propMC.particles_defaultHistories;
    matRad_cfg.dispInfo('Using default number of Histories per Bixel: %d\n',nCasePerBixel);
end

if nargin < 6
    calcDoseDirect = false;
end

if nargin < 7
    openStack = false;
end

if isfield(pln,'propMC') && isfield(pln.propMC,'outputVariance')
    matRad_cfg.dispWarning('Variance scoring for TOPAS not yet supported.');
end

if ~calcDoseDirect
    %     matRad_cfg.dispError('matRad so far only supports direct dose calculation for TOPAS!\n');
    matRad_cfg.dispWarning('You have selected TOPAS dij calculation, this may take a while ^^');
end

if ~isfield(pln.propStf,'useRangeShifter') 
    pln.propStf.useRangeShifter = false;
end

%if pln.propStf.useRangeShifter
%    matRad_cfg.dispError('matRad''s TOPAS interface does not support range shifters yet!\n');
%end

env = matRad_getEnvironment();

%% Initialize dose Grid as usual
matRad_calcDoseInit;

% fill bixels, rays and beams in case of dij calculation
if ~calcDoseDirect
    counter = 1;
    for f = 1:dij.numOfBeams
        for r = 1:stf(f).numOfRays
            for b = 1:stf(f).numOfBixelsPerRay(r)
                dij.bixelNum(counter) = b;
                dij.rayNum(counter)   = r;
                dij.beamNum(counter)  = f;
                counter = counter + 1;
            end
        end
    end
end

% for TOPAS we explicitly downsample the ct to the dose grid (might not
% be necessary in future versions with separated grids)


for s = 1:ct.numOfCtScen
    cubeHUresampled{s} =  matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y',  dij.ctGrid.z,ct.cubeHU{s}, ...
                                dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'linear');
    cubeResampled{s} =  matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y',  dij.ctGrid.z,ct.cube{s}, ...
                                dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'linear');                        
end

%Allocate temporary resampled CT
ctR = ct;
ctR.cube = cell(1);
ctR.cubeHU = cell(1);
ctR.numOfCtScen = 1;
ctR.resolution = dij.doseGrid.resolution;
ctR.cubeDim = dij.doseGrid.dimensions;
ctR.x = dij.doseGrid.x;
ctR.y = dij.doseGrid.y;
ctR.z = dij.doseGrid.z;

%% sending data to topas

load([pln.radiationMode,'_',pln.machine]);
topasConfig = MatRad_TopasConfig();
% Create Base Data
if isstruct(machine.data(i).Z)
    for i = 1:length(machine.data)
        machine.data(i).Z = machine.data(i).Z.profileORG;
    end
end
topasBaseData = MatRad_TopasBaseData(machine,stf);%,TopasConfig);

topasConfig.numHistories = nCasePerBixel;
if openStack
    topasConfig.workingDir = [topasConfig.thisFolder filesep 'MCrun' filesep [pln.machine,'_',pln.radiationMode,'_']];
    topasConfig.workingDir = [topasConfig.workingDir num2str(length(dir([topasConfig.workingDir,'*'])) + 1) filesep];
end
% topasConfig.numOfRuns = matRad_cfg.propMC.topas_defaultNumBatches;

%Collect weights
if calcDoseDirect
    w = zeros(sum([stf(:).totalNumOfBixels]),ct.numOfCtScen);
    counter = 1;
    for i = 1:length(stf)
        for j = 1:stf(i).numOfRays
            rayBix = stf(i).numOfBixelsPerRay(j);
            w(counter:counter+rayBix-1,:) = stf(i).ray(j).weight;
            counter = counter + rayBix;
        end
    end
end

currDir = cd;

for shiftScen = 1:pln.multScen.totNumShiftScen
    
    % manipulate isocenter
    for k = 1:length(stf)
        stf(k).isoCenter = stf(k).isoCenter + pln.multScen.isoShift(shiftScen,:);
    end    
    
    for ctScen = 1:pln.multScen.numOfCtScen
        for rangeShiftScen = 1:pln.multScen.totNumRangeScen
            if pln.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)
                
                %Overwrite CT (TEMPORARY - we should use 4D calculation in
                %TOPAS here)
                ctR.cubeHU = cubeHUresampled(ctScen);
                ctR.cube = cubeResampled(ctScen);
                %Delete previous topas files
                delete([topasConfig.workingDir,'score*'])
                delete([topasConfig.workingDir,'matRad_plan*'])

                if calcDoseDirect
                    topasConfig.writeAllFiles(ctR,pln,stf,topasBaseData,w(:,ctScen));
                else
                    topasConfig.writeAllFiles(ctR,pln,stf,topasBaseData);
                end
                
                % Run simulation for current scenario
                cd(topasConfig.workingDir);
                
                if openStack
                    save('dij.mat','dij')
                    save('weights.mat','w')
                    matRad_cfg.dispInfo('TOPAS simulation skipped for external calculation\n');
                else
                    
                    for beamIx = 1:numel(stf)
                        
                        for runIx = 1:topasConfig.numOfRuns
                            fname = sprintf('%s_field%d_run%d',topasConfig.label,beamIx,runIx);
                            topasCall = sprintf('%s %s.txt > %s.out > %s.err',topasConfig.topasExecCommand,fname,fname,fname);
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
                                fin = cellfun(@(f) exist(f,'file'),finishedFiles);
                                runsFinished = all(fin);
                            end
                        end
                        
                    end
                    
                    cd(currDir);
                    
                    %% Simulation finished - read out volume scorers from topas simulation
                    if calcDoseDirect
                        topasCubes = matRad_readTopasData(topasConfig.workingDir);
                    else
                        topasCubes = matRad_readTopasData(topasConfig.workingDir,dij);
                    end
                    
                    fnames = fieldnames(topasCubes);
                    dij.MC_tallies = fnames;
                    
                    if calcDoseDirect
                        if ~isfield(topasCubes,'RBE')
                            for f = 1:numel(fnames)
                                dij.(fnames{f}){ctScen,1} = sum(w(:,ctScen))*reshape(topasCubes.(fnames{f}),[],1);
                            end
                        else
                            for d = 1:length(stf)
                                dij.physicalDose{ctScen,1}(:,d)    = sum(w)*reshape(topasCubes.(['physicalDose_beam',num2str(d)]),[],1);
                                dij.alpha{ctScen,1}(:,d)           = reshape(topasCubes.(['alpha_beam',num2str(d)]),[],1);
                                dij.beta{ctScen,1}(:,d)            = reshape(topasCubes.(['beta_beam',num2str(d)]),[],1);
                                %                             dij.RBE{1}(:,d)             = reshape(topasCubes.(['RBE_beam',num2str(d)]),[],1);
                                
                                [dij.ax,dij.bx] = matRad_getPhotonLQMParameters(cst,prod(ct.cubeDim),1);
                                %                             dij.ax = full(reshape(ax,ct.cubeDim));
                                %                             dij.bx = full(reshape(bx,ct.cubeDim));
                                dij.abx(dij.bx>0) = dij.ax(dij.bx>0)./dij.bx(dij.bx>0);
                                
                                dij.mAlphaDose{ctScen,1}(:,d)      = dij.physicalDose{ctScen,1}(:,d) .* dij.alpha{ctScen,1}(:,d);
                                dij.mSqrtBetaDose{ctScen,1}(:,d)   = sqrt(dij.physicalDose{ctScen,1}(:,d)) .* dij.beta{ctScen,1}(:,d);
                            end
                        end
                    else
                        for f = 1:numel(fnames)
                            for d = 1:stf(f).totalNumOfBixels
                                dij.physicalDose{1}(:,d) = reshape(topasCubes.(fnames{f}){d},[],1);
                            end
                        end
                    end
                end
            end
        end
    end
end
    
    % manipulate isocenter back
    for k = 1:length(stf)
        stf(k).isoCenter = stf(k).isoCenter - pln.multScen.isoShift(shiftScen,:);
    end
end
