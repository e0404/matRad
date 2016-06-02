function ct = matRad_calcWaterEqD(ct)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to calculate the equivalent densities from a 
% dicom ct that originally uses intensity values
%
% call
%   ct = matRad_calcWaterEqD(ct)
%
% input
%   ct: unprocessed dicom ct data which are stored as intensity values (IV)
%
%                      HU = IV * slope + intercept
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
ctHU = double(ct.cube{1}) * double(ct.dicomInfo.RescaleSlope) + double(ct.dicomInfo.RescaleIntercept);

%% conversion from HU to water equivalent density

% load hlut

% directory with look up table files
hlutDir = fullfile(fileparts(mfilename('fullpath')),'hlutLibrary',filesep);

% if possible -> file standard out of dicom tags
try
    manufacturer = ct.dicomInfo.Manufacturer;
    model = ct.dicomInfo.ManufacturerModelName;
    convKernel = ct.dicomInfo.ConvolutionKernel;
    
    hlutFileName = strcat(manufacturer, '-', model, '-ConvolutionKernel-',...
        convKernel, '.hlut');
    
    % check whether fileNames used '-' or '_' instead of blanks
    hlutFileCell{1} = hlutFileName;
    hlutFileCell{2} = regexprep(hlutFileName,' ','-');
    hlutFileCell{3} = regexprep(hlutFileName,' ','_');
    
    % add pathname
    hlutFileCell = strcat(hlutDir,hlutFileCell);

    for i = 1:3
        existIx(i) = exist(hlutFileCell{i}, 'file') == 2;
    end

    if sum(existIx) == 0
        warnText = {['Could not find HLUT ' hlutFileName ' in hlutLibrary folder.' ...
            ' matRad default HLUT loaded']};
        warndlg(warnText,'Could not load HLUT');
        warning('matRad default HLUT loaded');
        % load default HLUT
        hlutFileName = strcat(hlutDir,'matRad_default.hlut');
    else
        hlutFileName = hlutFileCell{existIx};
    end

catch
    warnText = {['Could not construct hlut file name from DICOM tags.' ...
        ' matRad default HLUT loaded']};
    warndlg(warnText,'Could not load HLUT');
    warning('matRad default HLUT loaded');
       
    hlutFileName = strcat(hlutDir,'matRad_default.hlut');

end

hlutFile = fopen(hlutFileName,'r');
hlut = cell2mat(textscan(hlutFile,'%f %f','CollectOutput',1,'commentStyle','#'));
fclose(hlutFile);

% Manual adjustments if ct data is corrupt. If some values are out of range
% of the LUT, then these values are adjusted.
if max(ctHU(:)) > max(hlut(:,1))
    warning('projecting out of range HU values');
    ctHU(ctHU > max(hlut(:,1))) = max(hlut(:,1));
end
if min(ctHU(:)) < min(hlut(:,1))
    warning('projecting out of range HU values');
    ctHU(ctHU < min(hlut(:,1))) = min(hlut(:,1));
end

% interpolate HU to relative electron density based on lookup table
ct.cube{1} = interp1(hlut(:,1),hlut(:,2),double(ctHU));

% save hlut
ct.hlut = hlut;

end
