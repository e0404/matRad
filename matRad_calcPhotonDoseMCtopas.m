function dij = matRad_calcPhotonDoseMCtopas(ct,stf,pln,cst,calcDoseDirect)
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

matRad_cfg = MatRad_Config.instance();

if nargin < 5
    calcDoseDirect = false;
end

% Set parameters for full Dij calculation
if ~calcDoseDirect
    pln.propMC.scorer.calcDij = true;
    pln.propMC.numOfRuns = 1;
    
    % Load class variables in pln
    % for calcDoseDirect, this is already done in superior function
    if ~isa(pln.propMC,'MatRad_TopasConfig')
        pln = matRad_cfg.getDefaultClass(pln,'propMC','MatRad_TopasConfig');
    end
else
    if ~isa(pln.propMC,'MatRad_TopasConfig')
        matRad_cfg.dispError('Run calcParticleDoseMCtopas through calcDoseDirectMC');
    end
end

% load default parameters for doseCalc in case they haven't been set yet
pln = matRad_cfg.getDefaultProperties(pln,'propDoseCalc');

% set nested folder structure if external calculation is turned on (this will put new simulations in subfolders)
if pln.propMC.externalCalculation
    pln.propMC.workingDir = [pln.propMC.thisFolder filesep 'MCrun' filesep];
    pln.propMC.workingDir = [pln.propMC.workingDir pln.radiationMode,'_',pln.machine,'_',datestr(now, 'dd-mm-yy')];
end

%% Initialize dose Grid as usual
matRad_calcDoseInit;

% for TOPAS we explicitly downsample the ct to the dose grid (might not be necessary in future versions with separated grids)
[ctR,~,~] = matRad_resampleCTtoGrid(ct,cst,pln,stf);

% overwrite CT grid in dij in case of modulation.
if isfield(ctR,'ctGrid')
    dij.ctGrid = ctR.ctGrid;
end

%% sending data to topas

load([pln.radiationMode,'_',pln.machine]);

%Collect weights
if calcDoseDirect
    w = zeros(sum([stf(:).totalNumOfBixels]),1);
    ct = 1;
    for i = 1:length(stf)
        for j = 1:stf(i).numOfRays
            rayBix = stf(i).numOfBixelsPerRay(j);
            w(ct:ct+rayBix-1) = stf(i).ray(j).weight;
            ct = ct + rayBix;
        end
    end
else
    w = ones(sum([stf(:).totalNumOfBixels]),1);
end

currDir = cd;

for shiftScen = 1:pln.multScen.totNumShiftScen
    
    % manipulate isocenter
    for k = 1:length(stf)
        stf(k).isoCenter = stf(k).isoCenter + pln.multScen.isoShift(shiftScen,:);
    end
    
    % Delete previous topas files so there is no mix-up
    files = dir([pln.propMC.workingDir,'*']);
    files = {files(~[files.isdir]).name};
    fclose('all');
    for i = 1:length(files)
        delete([pln.propMC.workingDir,files{i}])
    end
    
    for ctScen = 1:pln.multScen.numOfCtScen
        for rangeShiftScen = 1:pln.multScen.totNumRangeScen
            if pln.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)
                
                % Save ctScen and rangeShiftScen for file constructor
                if ct.numOfCtScen > 1
                    ctR.currCtScen = ctScen;
                    ctR.currRangeShiftScen = rangeShiftScen;
                end
                
                % actually write TOPAS files
                if calcDoseDirect
                    pln.propMC.writeAllFiles(ctR,cst,pln,stf,machine,w);
                else
                    pln.propMC.writeAllFiles(ctR,cst,pln,stf,machine);
                end
                
            end
        end
    end
    
    % change director back to original directory
    cd(pln.propMC.workingDir);
    
    % Skip local calculation and data readout with this parameter. All necessary parameters to read the data back in
    % later are stored in the MCparam file that is stored in the folder. The folder is generated in the working
    % directory and the matRad_plan*.txt file can be manually called with TOPAS.
    if pln.propMC.externalCalculation
        matRad_cfg.dispInfo(['TOPAS simulation skipped for external calculation\nFiles have been written to: "',replace(pln.propMC.workingDir,'\','\\'),'"']);
    else
        for ctScen = 1:ct.numOfCtScen
            for beamIx = 1:numel(stf)
                for runIx = 1:pln.propMC.numOfRuns
                    if ct.numOfCtScen > 1
                        fname = sprintf('%s_field%d_ct%d_run%d',pln.propMC.label,beamIx,ctScen,runIx);
                    else
                        fname = sprintf('%s_field%d_run%d',pln.propMC.label,beamIx,runIx);
                    end
                    
                    if isprop(pln.propMC,'verbosity') && strcmp(pln.propMC.verbosity,'full')
                        topasCall = sprintf('%s %s.txt',pln.propMC.topasExecCommand,fname);
                    else
                        topasCall = sprintf('%s %s.txt > %s.out 2> %s.log',pln.propMC.topasExecCommand,fname,fname,fname);
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
                            matRad_cfg.dispError('TOPAS segmentation fault: could also be caused from an outdated TOPAS version or Linux distribution');
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
    end
    
    % revert back to original directory
    cd(currDir);
    
end

%% Simulation(s) finished - read out volume scorers from topas simulation
% Skip readout if external files were generated
if ~pln.propMC.externalCalculation
    dij = pln.propMC.readFiles(pln.propMC.workingDir);
    
    % Order fields for easier comparison between different dijs
    dij = orderfields(dij);
else
    dij = [];
end

% manipulate isocenter back
for k = 1:length(stf)
    stf(k).isoCenter = stf(k).isoCenter - pln.multScen.isoShift(shiftScen,:);
end
end
