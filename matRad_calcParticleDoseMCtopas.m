function topasCubes = matRad_calcParticleDoseMCtopas(ct,stf,pln,cst,calcDoseDirect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad TOPAS Monte Carlo proton dose calculation wrapper
%   This calls a TOPAS installation (not included in matRad due to
%   licensing model of TOPAS) for MC simulation
%
% call
%   dij = matRad_calcParticleDoseMCtopas(ct,stf,pln,cst,calcDoseDirect)
%
% input
%   ct:                         matRad ct struct
%   stf:                        matRad steering information struct
%   pln:                        matRad plan meta information struct
%   cst:                        matRad cst struct
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

% Instance of MatRad_Config class
matRad_cfg = MatRad_Config.instance();

% handle inputs
if nargin < 5
    calcDoseDirect = false;
end

if ~isfield(pln.propStf,'useRangeShifter')
    pln.propStf.useRangeShifter = false;
end

% Set parameters for full Dij calculation
if ~calcDoseDirect
    pln.propMC.scorer.calcDij = true;
    pln.propMC.numOfRuns = 1;

    % Load class variables in pln
    pln = matRad_cfg.getDefaultClass(pln,'propMC','MatRad_TopasConfig');

    if pln.propMC.numHistories  < 1e10
        matRad_cfg.dispWarning('Selected TOPAS dij calculation with fewer than normal histories (default 1e10), make sure you want to continue.');
    else
        matRad_cfg.dispWarning('You have selected TOPAS dij calculation, this may take a while ^^');
    end
else
    if ~isa(pln.propMC,'MatRad_TopasConfig')
        matRad_cfg.dispError('Run calcParticleDoseMCtopas through calcDoseDirectMC');
    end
end

% load default parameters for doseCalc in case they haven't been set yet
pln = matRad_cfg.getDefaultProperties(pln,{'propDoseCalc'});

% set nested folder structure if external calculation is turned on (this will put new simulations in subfolders)
if pln.propMC.externalCalculation
    pln.propMC.workingDir = [pln.propMC.thisFolder filesep 'MCrun' filesep];
    if isfield(pln,'patientID')
        pln.propMC.workingDir = [pln.propMC.workingDir pln.radiationMode filesep pln.patientID filesep pln.patientID '_'];
    end
    pln.propMC.workingDir = [pln.propMC.workingDir pln.radiationMode,'_',pln.machine,'_',datestr(now, 'dd-mm-yy')];
    if isfield(ct,'sampleIdx')
        pln.propMC.workingDir = [pln.propMC.workingDir '_' num2str(ct.sampleIdx,'%02.f') filesep];
    end
end

%% Initialize dose grid and dij

% load calcDoseInit as usual
matRad_calcDoseInit;

% for TOPAS we explicitly downsample the ct to the dose grid (might not be necessary in future versions with separated grids)
[ctR,~,~] = pln.propMC.resampleGrid(ct,cst,pln,stf);

% overwrite CT grid in dij in case of modulation.
if isfield(ctR,'ctGrid')
    dij.ctGrid = ctR.ctGrid;
end

%% sending data to topas

% Load and create TOPAS Base Data
load([pln.radiationMode,'_',pln.machine],'machine');

% Collect given weights
if calcDoseDirect
    w = zeros(sum([stf(:).totalNumOfBixels]),ctR.numOfCtScen);
    counter = 1;
    for i = 1:length(stf)
        for j = 1:stf(i).numOfRays
            rayBix = stf(i).numOfBixelsPerRay(j);
            w(counter:counter+rayBix-1,:) = stf(i).ray(j).weight;
            counter = counter + rayBix;
        end
    end
end

% Get photon parameters for RBExD calculation
if isfield(pln,'bioParam') && strcmp(pln.bioParam.quantityOpt,'RBExD')
    pln.propMC.scorer.RBE = true;
    [dij.ax,dij.bx] = matRad_getPhotonLQMParameters(cst,dij.doseGrid.numOfVoxels,1,VdoseGrid);
    dij.abx(dij.bx>0) = dij.ax(dij.bx>0)./dij.bx(dij.bx>0);
end

% save current directory to revert back to later
currDir = cd;

for shiftScen = 1:pln.multScen.totNumShiftScen

    % manipulate isocenter
    for k = 1:length(stf)
        stf(k).isoCenter = stf(k).isoCenter + pln.multScen.isoShift(shiftScen,:);
    end

    % Run simulations for each scenario
    for ctScen = 1:pln.multScen.numOfCtScen
        for rangeShiftScen = 1:pln.multScen.totNumRangeScen
            if pln.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)

                % Delete previous topas files so there is no mix-up
                files = dir([pln.propMC.workingDir,'*']);
                files = {files(~[files.isdir]).name};
                fclose('all');
                for i = 1:length(files)
                    delete([pln.propMC.workingDir,files{i}])
                end

                % actually write TOPAS files
                if calcDoseDirect
                    pln.propMC.writeAllFiles(ctR,cst,pln,stf,machine,w(:,ctScen));
                else
                    pln.propMC.writeAllFiles(ctR,cst,pln,stf,machine);
                end

                % change director back to original directory
                cd(pln.propMC.workingDir);

                % save dij and weights, they are needed for later reading the data back in
                if pln.propMC.externalCalculation
                    matRad_cfg.dispInfo('TOPAS simulation skipped for external calculation\n');
                else
                    for beamIx = 1:numel(stf)
                        for runIx = 1:pln.propMC.numOfRuns
                            fname = sprintf('%s_field%d_run%d',pln.propMC.label,beamIx,runIx);
                            if isfield(pln.propMC,'verbosity') && strcmp(pln.propMC.verbosity,'full')
                                topasCall = sprintf('%s %s.txt',pln.propMC.topasExecCommand,fname);
                            else
                                topasCall = sprintf('%s %s.txt > %s.out > %s.log',pln.propMC.topasExecCommand,fname,fname,fname);
                            end

                            % initiate parallel runs and delete previous files
                            if pln.propMC.parallelRuns
                                finishedFiles{runIx} = sprintf('%s.finished',fname);
                                topasCall = [topasCall '; touch ' finishedFiles{runIx} ' &'];
                            end

                            % Actual simulation happening here
                            matRad_cfg.dispInfo('Calling TOPAS: %s\n',topasCall);
                            [status,cmdout] = system(topasCall,'-echo');

                            % Process TOPAS output and potential errors
                            cout = splitlines(string(cmdout));
                            if status == 0
                                matRad_cfg.dispInfo('TOPAS simulation completed succesfully\n');
                            else
                                if status == 139
                                    matRad_cfg.dispError('TOPAS segmentation fault: might be caused from an outdated TOPAS version or Linux distribution');
                                else
                                    matRad_cfg.dispError('TOPAS simulation exited with error code %d\n "%s"',status,cout(2:end-1));
                                end
                            end
                        end

                        % wait for parallel runs to finish and process
                        if pln.propMC.parallelRuns
                            runsFinished = false;
                            pause('on');
                            while ~runsFinished
                                pause(1);
                                fin = cellfun(@(f) exist(f,'file'),finishedFiles);
                                runsFinished = all(fin);
                            end
                            % Delete marker files
                            delete(finishedFiles{:});
                        end

                    end
                end

                % revert back to original directory
                cd(currDir);

            end
        end
    end
end

%% Simulation(s) finished - read out volume scorers from topas simulation
% Skip readout if external files were generated
if ~pln.propMC.externalCalculation
    topasCubes = pln.propMC.readFiles(pln.propMC.workingDir);
else
    topasCubes = [];
end

% manipulate isocenter back
for k = 1:length(stf)
    stf(k).isoCenter = stf(k).isoCenter - pln.multScen.isoShift(shiftScen,:);
end
end
