function [resultGUI] = matRad_importDicomRTDose(ct, rtDoseFiles, pln)
% matRad function to import dicom RTDOSE data
% 
% call
%   resultGUI = matRad_importDicomRTDose(ct, rtDoseFiles)
%   resultGUI = matRad_importDicomRTDose(ct, rtDoseFiles, pln)
%
% input
%   ct:             ct imported by the matRad_importDicomCt function
%   rtDoseFiles:   	cell array of RTDOSE Dicom files
%   pln:            (optional) matRad pln struct
%
% output
%   resultGUI:      matRad resultGUI struct with different beams. Note that
%                   the summation (called plan) of the beams is named 
%                   without subscripts, e.g. physical_Dose.
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
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

%% import and interpolate dose files
% number of dosefiles
numDoseFiles = size(rtDoseFiles,1);

for i = 1 : numDoseFiles
    currDose = rtDoseFiles(i,:);
    itemName = strcat('Item_',num2str(i));
    dose.(itemName) = matRad_interpDicomDoseCube( ct, currDose);
end

%% put dose information and dose meta information to resultGUI
countBeamNumberPhysDose = 1;
countBeamNumberRBExDose = 1;
countBeamNumberOther = 1;
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
    
    if strncmpi(doseTypeHelper,'PHYSICAL',6)
        doseTypeHelper = 'physicalDose';
    elseif strncmpi(doseTypeHelper,'EFFECTIVE',6)
        doseTypeHelper = 'RBExDose';
    end
    
    %If given as plan and not per fraction
    if strcmpi(doseSumHelper,'PLAN') || strcmpi(doseSumHelper,'BEAM')
        if exist('pln','var') 
            dose.(itemName).cube = dose.(itemName).cube / pln.numOfFractions;
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
    
    resultGUI.(resultName) = dose.(itemName).cube;
    resultGUI.doseMetaInfo.(resultName) = dose.(itemName).dicomInfo;
    
end
% save timeStamp
resultGUI.doseMetaInfo.timeStamp = datestr(clock);
end
