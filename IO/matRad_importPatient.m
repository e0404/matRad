function [ct,cst] = matRad_importPatient(ctFile,maskFiles,HLUT)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad patient import from binary files (CT and masks)
% 
% call
%   [ct,cst] = matRad_importPatient(cubeFile,maskFiles)
%
% input
%   ctFile:     path to CT file. If HLUT is not set, values are interpreted
%               as HU and interpolated to ED. 
%   maskFiles:  cell array with filenames to the masks
%               if maskFiels contains a folder, all contained and
%               recognized data files are treated as masks
%   HLUT:       optional HLUT, (n,2) array. if set to 'default', we will
%               use a default HLUT
% output
%   ct          ct struct for use with matlab
%   cst         cst struct for use with matlab
%
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


[cube, metadata] = matRad_readCube(ctFile);

ct = struct();
cst = cell(0,6);
%Create dummy cst?

if nargin < 3
    ct.cube{1} = cube;
else 
    if isequal(HLUT,'default')
            HLUT = [-1024.0 0.00324; ...
            200.0   1.20000; ...
            449.0   1.20000; ...
            2000.0  2.49066; ...	 	
            2048    2.53060; ...
            3071    2.53060];
    end
    ct.cube{1} = interp1(HLUT(:,1),HLUT(:,2),double(cube));
    ct.hlut = HLUT;
end

ct.cubeDim = metadata.cubeDim;
ct.resolution.x = metadata.resolution(1);
ct.resolution.y = metadata.resolution(2);
ct.resolution.z = metadata.resolution(3);
    
ct.numOfCtScen = 1;

maskId = 1;
for maskFile=maskFiles(:)'
    if exist(maskFile{1},'dir')
        contents = dir(maskFile{1});
        for s=1:numel(contents)
            if(~contents(s).isdir)
                [mask, maskMeta] = matRad_readCube(fullfile(maskFile{1},contents(s).name));
                cstLine = importMaskToCstLine(maskId,mask,maskMeta);
                cst = [cst; cstLine]; 
                maskId = maskId + 1;
            end
        end
    elseif exist(maskFile{1},'file')
        [mask,maskMeta] = matRad_readCube(maskFile{1});  
        cstLine = importMaskToCstLine(maskId,mask,maskMeta);
        cst = [cst; cstLine];           
        maskId = maskId + 1;
    else
        disp(['Ignored file/dir ' maskFile '!']);
    end
end

end

function cstLine = importMaskToCstLine(maskId,mask,maskMeta)  
    cstLine = cell(1,6);
    cstLine{1} = maskId - 1;
    cstLine{2} = maskMeta.name;
    cstLine{3} = 'OAR';
    cstLine{4}{1} = find(mask > 0);
    cstLine{5}.Priority = maskId;
    cstLine{5}.alphaX = 0.1;
    cstLine{5}.betaX = 0.05;
    cstLine{5}.Visible = 1;
end


