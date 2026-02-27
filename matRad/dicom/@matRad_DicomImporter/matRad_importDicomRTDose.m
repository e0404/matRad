function obj = matRad_importDicomRTDose(obj)
% matRad function to import dicom RTDOSE data
% 
% In your object, there must be properties that contain:
%   - ct imported by the matRad_importDicomCt function;
%   - cell array of RTDose DICOM files.
% Optional:
%   - matRad pln structure.
%
% Output - matRad resultGUI structure with different beams.
% Note that the summation (called plan) of the beams is named without 
% subscripts, e.g. physical_Dose.
%
% call
%   obj = matRad_importDicomRTDose(obj)
%
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

%% import and interpolate dose files
% number of dosefiles
numDoseFiles = size(obj.importFiles.rtdose,1);

for i = 1 : numDoseFiles
    obj.importRTDose.currDose = obj.importFiles.rtdose(i,:);
    itemName = strcat('Item_',num2str(i));
    obj = matRad_interpDicomDoseCube(obj);
    dose.(itemName) = obj.importRTDose.dose;
end

%% put dose information and dose meta information to resultGUI
countBeamNumberPhysDose = 1;
countBeamNumberRBExDose = 1;
countBeamNumberOther = 1;

obj.resultGUI = struct();

for i = 1 : numDoseFiles
    itemName = strcat('Item_',num2str(i));
    doseTypeHelper      = dose.(itemName).dicomInfo.DoseType;
    doseSumHelper       = dose.(itemName).dicomInfo.DoseSummationType;
    
    %Field is not always existing
    if isfield(dose.(itemName).dicomInfo,'InstanceNumber')
        doseInstanceHelper  = dose.(itemName).dicomInfo.InstanceNumber;
    else
        doseInstanceHelper = [];
    end

    doseComment = '';
    if isfield(dose.(itemName).dicomInfo,'DoseComment')
        doseComment = dose.(itemName).dicomInfo.DoseComment;
    end
            
    if strncmpi(doseTypeHelper,'PHYSICAL',6)
        if any(strcmp(doseComment,{'physicalDose','LET'}))
            doseTypeHelper = doseComment;
        else
            doseTypeHelper = 'physicalDose';
        end
    elseif strncmpi(doseTypeHelper,'EFFECTIVE',6)
        if any(strcmp(doseComment,{'RBExDose','BED','alpha','beta','effect','RBE'}))
            doseTypeHelper = doseComment;
        else
            doseTypeHelper = 'RBExDose';
        end
    end
    
    %If given as plan and not per fraction
    if strcmpi(doseSumHelper,'PLAN') || strcmpi(doseSumHelper,'BEAM')
        if ~isempty(obj.pln) 
            dose.(itemName).cube = dose.(itemName).cube / obj.pln.numOfFractions;
        else
            matRad_cfg.dispWarning('DICOM dose given as PLAN, but no pln struct available to compute fraction dose! Assuming 1 fraction!');
        end
    end
        
    if strncmpi(doseSumHelper,'BEAM',4)
        try
            beamNumber = dose(itemName).dicomInfo.ReferencedRTPlanSequence.ReferencedFractionGroupSequence.ReferencedBeamSequence.ReferencedBeamNumber;
        catch
            switch doseTypeHelper
                case 'physicalDose'
                    beamNumber = countBeamNumberPhysDose;
                    countBeamNumberPhysDose = countBeamNumberPhysDose +1;
                case 'RBExDose'
                    beamNumber = countBeamNumberRBExDose;
                    countBeamNumberRBExDose = countBeamNumberRBExDose + 1;
                otherwise
                    beamNumber = countBeamNumberOther;
                    countBeamNumberOther = countBeamNumberOther + 1;
            end
        end
        
        beamSuffix = ['_beam' num2str(beamNumber)];
    else
        beamSuffix = '';
    end
    
    if ~isempty(doseInstanceHelper)
        instanceSuffix = ['_' num2str(doseInstanceHelper)];
    else
        instanceSuffix = '';
    end        
    
    resultName = strcat(doseTypeHelper,instanceSuffix,beamSuffix);
        

    if isfield(obj.resultGUI,resultName)
        count = 1;
        addSuffix = ['_' num2str(count)];
        resultNameNew = [resultName addSuffix];
        while isfield(obj.resultGUI,resultNameNew)
            count = count + 1;
            addSuffix = ['_' num2str(count)];
            resultNameNew = [resultName addSuffix];
        end       
               
        matRad_cfg.dispWarning('Already imported dose ''%s'', naming the new dose ''%s'', manually organize resultGUI afterwards for duplicates!',resultName, resultNameNew);
        resultName = resultNameNew;
    end
    
    obj.resultGUI.(resultName) = dose.(itemName).cube;
    obj.resultGUI.doseMetaInfo.(resultName) = dose.(itemName).dicomInfo;
    
end
% save timeStamp
obj.resultGUI.doseMetaInfo.timeStamp = datestr(clock);
end
