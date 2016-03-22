function  [ Flag ]=  matRad_unitTestGenData(plnNew,sPatientName)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to create a reference unit test result
% 
% call
%   matRad_unitTestGenData(plnNew,sPatientName)
%
% input
%   plnNew:       matRad's pln struct 
%   sPatientName: String which defines the patient data e.g. 'BOXPHANTOM'
%
% output
%   Flag:         boolean if the data set generation was successfull
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

Flag = false;
plnNew.patientname = sPatientName;
%dirty
rootDir = pwd;
rootDir = rootDir(1:end-9);  % /unitTest has exactly nine chars
addpath(rootDir);
sParticle = '';
switch plnNew.radiationMode
    case 'photons'
        sParticle = 'X';
    case 'protons'
        sParticle = 'P';
    case 'carbon'
        sParticle = 'C';
end

%% check how many plans for this patient data exist 
sNamePrefix = ['unitTest_result_' sParticle '_'];
sNamePostfix = ['_' sPatientName '.mat'];
ExistingPlans = dir([pwd filesep sNamePrefix '*.mat']);
plnID = 1;

% loop over all existing unit test results
for i = 1:length(ExistingPlans)
    % find results which have the selected Patientname
    if strfind(ExistingPlans(i).name,sPatientName)>0
        % load plan and check if they are the same - if so ask for
        % replacement
        load(ExistingPlans(i).name)
        
        if ((pln.numOfFractions  == plnNew.numOfFractions)     && ...
            strcmp(pln.bioOptimization,plnNew.bioOptimization) && ...
            strcmp(pln.radiationMode,plnNew.radiationMode)     && ...
            pln.numOfBeams      == plnNew.numOfBeams           && ...
            isequal(pln.couchAngles,plnNew.couchAngles)        && ...
            isequal(pln.gantryAngles,plnNew.gantryAngles)      && ...
            pln.bixelWidth      == plnNew.bixelWidth           && ...
           
                choice = questdlg('This plan already exists. Would you like to replace it ?',...
                    'Replace Plan',...
                    'Yes-Replace','No-create new','Abort','Yes-Replace');

                switch choice 
                    case 'Yes-Replace'
                        break
                    case 'Abort'
                        return
                    case 'No-create new'
                end
                  
            
        end
        
        plnID = plnID + 1;  
    end
end

%% create result / either for test or for creation
load([sPatientName])
plnNew.isoCenter = matRad_getIsoCenter(cst,ct,0);
plnNew.numOfVoxels     = numel(ct.cube);
plnNew.voxelDimensions = size(ct.cube);

stf = matRad_generateStf(ct,cst,plnNew);

if strcmp(plnNew.radiationMode,'photons')
    dij = matRad_calcPhotonDose(ct,stf,plnNew,cst,0);
elseif strcmp(plnNew.radiationMode,'protons') || ...
        strcmp(plnNew.radiationMode,'carbon')
    dij = matRad_calcParticleDose(ct,stf,plnNew,cst,0);
end

resultGUI = matRad_fluenceOptimization(dij,cst,plnNew,0);

%% sequencing
if strcmp(plnNew.radiationMode,'photons') && (plnNew.runSequencing || plnNew.runDAO)
    %resultGUI = matRad_xiaLeafSequencing(resultGUI,stf,dij,5);
    resultGUI = matRad_engelLeafSequencing(resultGUI,stf,dij,5);
end

%% DAO
if strcmp(plnNew.radiationMode,'photons') && plnNew.runDAO
   resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,1);
   matRad_visApertureInfo(resultGUI.apertureInfo);
end

%% extract dij from iso center slice according to the sample spacing
plnNew.dijSpacing = 2;
resultGUI.dijIsoCenter = matRad_unitTestGetIsoCenterSlice(plnNew,dij,plnNew.dijSpacing);

%% save expected result
if plnID < 10
     plnID = ['00' num2str(plnID)];
elseif plnID <100
     plnID = ['0' num2str(plnID)];
else
     plnID = num2str(plnID);
end
pln = plnNew;
ExpResultGUI = resultGUI;
save([sNamePrefix plnID sNamePostfix],'pln', 'ExpResultGUI');

Flag = true;  

end



function dijIsoCenter = matRad_unitTestGetIsoCenterSlice(pln,dij,Spacing)

StartIdx = pln.voxelDimensions(1)*pln.voxelDimensions(2)*(pln.voxelDimensions(3)-1);
EndIdx   = pln.voxelDimensions(1)*pln.voxelDimensions(2)*(pln.voxelDimensions(3)); 
LinearIdx = StartIdx:Spacing:EndIdx;
dijIsoCenter = dij.physicalDose(LinearIdx,:);

end
