function ct = matRad_electronDensitiesToHU(ct)

% load hlut

% directory with look up table files
hlutDir = fullfile(fileparts(mfilename('fullpath')),'dicomImport',filesep,'hlutLibrary',filesep);

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

% interpolate rel. electron dens. to HU based on lookup table
if isequal(hlut(:,2),unique(hlut(:,2))) && isequal(hlut(:,1),unique(hlut(:,1)))
    ct.cubeHU{1} = interp1(hlut(:,2),hlut(:,1),ct.cube{1});
else
    error('Conversion from electron densities to HU values could not be done because HLUT is not bijective.');
end