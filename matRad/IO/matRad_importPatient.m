function [ct,cst] = matRad_importPatient(ctFile,maskFiles,hlutFilename)
% matRad patient import from binary files (CT and masks)
% 
% call
%   [ct,cst] = matRad_importPatient(cubeFile,maskFiles)
%   [ct,cst] = matRad_importPatient(cubeFile,maskFiles, hlutFilename)
%
% input
%   ctFile:     path to CT file. If HLUT is not set, values are interpreted
%               as HU and interpolated to ED. 
%   maskFiles:  cell array with filenames to the masks
%               if maskFiels contains a folder, all contained and
%               recognized data files are treated as masks
%   hlutFilname:(optional) HLUT, (n,2) array. if set to 'default', we will
%               use a default HLUT
% output
%   ct          ct struct for use with matlab
%   cst         cst struct for use with matlab
%
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


[cube, metadata] = matRad_readCube(ctFile);

ct = struct();
cst = cell(0,6);

ct.cubeHU{1} = cube;

if nargin == 3
    HLUT = matRad_readHLUT(hlutFilename);
    ct.cube{1} = interp1(HLUT(:,1),HLUT(:,2),double(cube));
    ct.hlut = HLUT;
end

ct.cubeDim = metadata.cubeDim;
ct.resolution.x = metadata.resolution(1);
ct.resolution.y = metadata.resolution(2);
ct.resolution.z = metadata.resolution(3);
    
ct.numOfCtScen = 1;

maskId = 1;
hGlobalWaitbar = waitbar(0,'Importing Segmentations');
set(findall(hGlobalWaitbar,'type','text'),'Interpreter','none');



for f=1:numel(maskFiles)
    maskFile = maskFiles{f};
    waitbar(f/numel(maskFiles),hGlobalWaitbar,['Importing Segmentations: ' maskFiles{f}]);
    if exist(maskFile,'dir')
        contents = dir(maskFile);
        hFolderWaitbar = waitbar(0,'Importing Folder');
        set(findall(hFolderWaitbar,'type','text'),'Interpreter','none');
        for s=1:numel(contents)
            waitbar(s/numel(contents),hFolderWaitbar,['Importing Folder: ' contents(s).name]);            
            if(~contents(s).isdir)
                [mask, maskMeta] = matRad_readCube(fullfile(maskFile,contents(s).name));
                cstLine = importMaskToCstLine(maskId,mask,maskMeta);
                cst = [cst; cstLine]; 
                maskId = maskId + 1;
            end            
        end
        delete(hFolderWaitbar);
    elseif exist(maskFile,'file')
        [mask,maskMeta] = matRad_readCube(maskFile);  
        cstLine = importMaskToCstLine(maskId,mask,maskMeta);
        cst = [cst; cstLine];           
        maskId = maskId + 1;
    else
        disp(['Ignored file/dir ' maskFile '!']);
    end
end

delete(hGlobalWaitbar);


%Assign default colors
colors = colorcube(size(cst,1));
for i = 1:size(cst,1)
    cst{i,5}.visibleColor = colors(i,:);
end

end

function cstLine = importMaskToCstLine(maskId,mask,maskMeta)  
    cstLine = cell(1,6);
    cstLine{1} = maskId - 1;
    cstLine{2} = maskMeta.name;
    cstLine{3} = tryToGetVoiTypeByName(maskMeta.name);
    cstLine{4}{1} = find(mask > 0);
    cstLine{5}.Priority = maskId;
    cstLine{5}.alphaX = 0.1;
    cstLine{5}.betaX = 0.05;
    cstLine{5}.Visible = 1;
end

function type = tryToGetVoiTypeByName(voiName)
    targetNames = {'target'; 'ptv'; 'ctv'};
    for n=1:numel(targetNames)
        found = strfind(lower(voiName),lower(targetNames{n}));
        if ~isempty(found)
            type = 'TARGET';
            return;
        end
    end
    type = 'OAR';
end


