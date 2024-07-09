function [cube, metadata] = matRad_readNifTI(filename)
% matRad NifTI reader
% 
% call
%   [cube, metadata] = matRad_readNifTI(filename)
%
% input
%   filename:   full path to .nii(.gz) file
%
% output
%   cube:       the read cube
%   metadata:   metadata from header information
%
% References
%   [1] http://teem.sourceforge.net/nrrd/format.html5
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2024 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%Instantiate matrad config
matRad_cfg = MatRad_Config.instance();

%Check for Image Processing toolbox
avail = license('test', 'image_toolbox');
if ~avail
    matRad_cfg.dispError('Image Processing Toolbox is required for reading NifTI files!\n');
else
    success = license('checkout', 'image_toolbox');
    if ~success
        matRad_cfg.dispError('Image Processing Toolbox license could not be checked out!\n');
    end
end

%Obtain NifTI info strcuture
try
    info = niftiinfo(filename);
catch ME
    matRad_cfg.dispError('Error reading NifTI file: %s\n', ME.message);
end

% Read the cube
try
    cube = niftiread(info);
catch ME
    matRad_cfg.dispError('Error reading NifTI file: %s\n', ME.message);
end

metadata = struct();
metadata.datatype = info.Datatype;
metadata.cubeDim = info.ImageSize;
metadata.axisPermutation = [1 2 3];
metadata.transform = info.Transform.T;
metadata.resolution = info.PixelDimensions;

