function [ resultGUI ] = matRad_importDicomRTDose(ct, rtDoseFiles, pln)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to import dicom RTDOSE data
% 
% call
%   resultGUI = matRad_importDicomRTDose(ct, rtDoseFiles)
%
% input
%   ct:             ct imported by the matRad_importDicomCt function
%   rtDoseFiles:   	cell array of RTDOSE Dicom files
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


%% import and interpolate dose files
% number of dosefiles
numDoseFiles = size(rtDoseFiles,1);

for i = 1 : numDoseFiles
    currDose = rtDoseFiles(i,:);
    itemName = strcat('Item_',num2str(i));
    dose.(itemName) = matRad_interpDicomDoseCube( ct, currDose);
end

%% put dose information and dose meta information to resultGUI
for i = 1 : numDoseFiles
    itemName = strcat('Item_',num2str(i));
    doseTypeHelper      = dose.(itemName).dicomInfo.DoseType;
    doseSumHelper       = dose.(itemName).dicomInfo.DoseSummationType;
    doseInstanceHelper  = dose.(itemName).dicomInfo.InstanceNumber;
    
    if strncmpi(doseTypeHelper,'PHYSICAL',6)
        doseTypeHelper = 'physicalDose_';
    elseif strncmpi(doseTypeHelper,'EFFECTIVE',6)
        doseTypeHelper = 'RBExDose_';
    end
    
    resultName = strcat(doseTypeHelper,doseSumHelper,'_',num2str(doseInstanceHelper));
    if strncmpi(resultName,'physicalDose_PLAN',17)
        resultName = 'physicalDose';
    elseif strncmpi(resultName,'RBExDose_PLAN',10)
        resultName = 'RBExDose';
    end
    
    % scale to fraction based dose
    if exist('pln','var')
        dose.(itemName).cube = dose.(itemName).cube / pln.numOfFractions;
    end
    
    resultGUI.(resultName) = dose.(itemName).cube;
    resultGUI.doseMetaInfo.(resultName) = dose.(itemName).dicomInfo;
    
end
% save timeStamp
resultGUI.doseMetaInfo.timeStamp = datestr(clock);
end