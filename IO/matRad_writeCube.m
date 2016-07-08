function matRad_writeCube(filepath,cube,datatype,metadata)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad wrapper for Cube export
% 
% call
%   matRad_writeCube(filepath,cube,meta)
%
% input
%   filepath:                   full output path. needs the right extension
%                               to choose the appropriate writer
%   cube:                       cube that is to be written
%   meta:                       meta-information in struct. 
%                               Necessary fieldnames are:
%                               - resolution: [x y z]
%                               Optional:
%                               - axisPermutation (matRad default [2 1 3])
%                               - coordinateSystem (matRad default 'LPS')
%                               - imageOrigin (as used in DICOM)
%                               - compress (true/false)
%                                 (default chosen by writer)
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

%% Sanity check
[filedir,filename,ext] = fileparts(filepath);
if ~exist(filedir,'dir')
    error(['Directory ' filedir ' does not exist!']);
end

%% Prepare Metadata
%if the field is not set, we assume standard matRad x-y swap
if ~isfield(metadata,'axisPermutation')
    metadata.axisPermutation = [2 1 3]; %Default Matlab axis permutation
end
%use the matrad coordinate system
if ~isfield(metadata,'coordinateSystem')
    metadata.coordinateSystem = 'LPS';  %Matlab coordinate system
end
%If there is no image origin set, center the image
if ~isfield(metadata,'imageOrigin')
    imageExtent = metadata.resolution .* size(cube);
    metadata.imageOrigin = zeros(1,numel(imageExtent)) - (imageExtent/2);
end
    
%% Choose writer
%So far we only have an nrrd writer
switch ext
    case '.nrrd'
        matRad_writeNRRD(filepath,cube,datatype,metadata);
    otherwise
        errordlg(['No writer found for extension "' ext '"']);
end


end

