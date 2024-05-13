function hlut = matRad_loadHLUT(ct, radiationMode)
% matRad function to load HLUT file based on the provided ct
%
% call
%   hlut = matRad_loadHLUT(ct, pln)
%
% input
%   ct:             ct with dicom information
%   radiationMode:  radiationMode as character array (e.g. 'photons') since matRad 3.
%                   Can also be a pln-struct for downwards compatibility
%
% output
%   hlut: lookup table
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team. 
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
  
% directory with look up table files
if ~isdeployed
    hlutDir = fullfile(fileparts(mfilename('fullpath')),'hlutLibrary',filesep);
else
    hlutDir = [];
end

%Old version - takes radiation mode from pln
if isstruct(radiationMode) && isfield(radiationMode,'radiationMode')
    particle = radiationMode.radiationMode;
elseif ischar(radiationMode)
    particle = radiationMode;
else
    matRad_cfg.dispError('Invalid radiation mode!');
end


% if possible -> file standard out of dicom tags
try
    
    hlutFileName = '';
    manufacturer = ct.dicomInfo.Manufacturer;
    model        = ct.dicomInfo.ManufacturerModelName;
    if isfield(ct.dicomInfo,'ConvolutionKernel')
        convKernel = ct.dicomInfo.ConvolutionKernel;
    else
        convKernel = '0';
    end
    hlutFileName = strcat(manufacturer, '-', model, '-ConvolutionKernel-',...
        convKernel, '_', particle, '.hlut');
    
    % check whether fileNames used '-' or '_' instead of blanks
    hlutFileCell{1} = hlutFileName;
    hlutFileCell{2} = regexprep(hlutFileName,' ','-');
    hlutFileCell{3} = regexprep(hlutFileName,' ','_');
    
    % add pathname
    hlutFileCell = strcat(hlutDir,hlutFileCell);

    % check if files exist
    existIx = cellfun(@(x) exist(x,'file') == 2,hlutFileCell);
    
    if sum(existIx) == 0
        produce an error to enter catch statment below :)
    else
        hlutFileName = hlutFileCell{existIx};
    end

catch
    
    warnText = ['Could not find HLUT ' hlutFileName ' in hlutLibrary folder.' ...
                ' matRad default HLUT loaded'];
    warning('off','backtrace')
    matRad_cfg.dispWarning(warnText);
    
    % load default HLUT
    hlutFileName = strcat(hlutDir,'matRad_default.hlut');

end

hlut = matRad_readHLUT(hlutFileName);

end
