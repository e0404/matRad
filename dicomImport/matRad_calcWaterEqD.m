function ct = matRad_calcWaterEqD(ct)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to calculate the equivalent densities from a 
% dicom ct that originally uses intensity values
%
% call
%   ct = matRad_calcWaterEqD(ct, slope, intercept)
%
% input
%   ct:                 unprocessed dicom ct data which are stored as
%                       intensity values (IV)
%
%                      HU = IV * slope + intercept
%
% base data
%   HU2waterEqT.mat:    look up table to convert from HU to relative 
%                       electron densities. Note that this is just example
%                       data. If you want to make precise computaitons for
%                       your scanner you need to replace this lookup table
%
% output
%   ct:                 ct struct with cube with relative _electron_ densities 
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


%% conversion from IV to HU
ctHU = double(ct.cube) * double(ct.dicomInfo.RescaleSlope) + double(ct.dicomInfo.RescaleIntercept);

%% conversion from HU to water equivalent density
%  load hlut
%  file standard
%  out of dicom tags
manufacturer = ct.dicomInfo_org.Manufacturer;
model = ct.dicomInfo_org.ManufacturerModelName;
convKernel = ct.dicomInfo_org.ConvolutionKernel;

clear hlutFileName hlutFile hlut
% convention of fileName
hlutFileName = strcat(manufacturer, '-', model, '-ConvolutionKernel-',...
    convKernel, '.hlut');

% check whether fileNames used '-' or '_' instead of blanks
hlutFileCell{1} = hlutFileName;
hlutFileCell{2} = regexprep(hlutFileName,' ','-');
hlutFileCell{3} = regexprep(hlutFileName,' ','_');

for i = 1:3
    if exist(hlutFileCell{i}, 'file') == 2
        existIx(i) = 1;
    else
        existIx(i) = 0;
    end
end

if sum(existIx) ~= 1
    warnText={'No or more than one proper HLUT corresponding to the DICOM tags loaded';...
        'Please provide .hlut file in hlutLibrary folder'};
    warndlg(warnText,'Could not load HLUT');
    % load default HLUT
    hlutFileName = 'matRad_default.hlut';
    hlutFile = fopen(hlutFileName,'r');
    hlut = cell2mat(textscan(hlutFile,'%f %f','CollectOutput',1,'commentStyle','#'));
    fclose(hlutFile);
    warning('matRad default HLUT loaded');
else
    hlutFileName = hlutFileCell{find(existIx)};
    hlutFile = fopen(hlutFileName,'r');
    hlut = cell2mat(textscan(hlutFile,'%f %f','CollectOutput',1,'commentStyle','#'));
    fclose(hlutFile);
end

if max(ctHU(:)) > max(hlut(:,1))
    warning('projecting out of range HU values');
    ctHU(ctHU>max(hlut(:,1))) = max(hlut(:,1));
end
if min(ctHU(:)) < min(hlut(:,1))
    warning('projecting out of range HU values');
    ctHU(ctHU<min(hlut(:,1))) = min(hlut(:,1));
end

% interpolate HU to relative electron density based on lookup table
ct.cube = interp1(hlut(:,1),hlut(:,2),double(ctHU));

% save hlut
ct.hlut = hlut;

end
