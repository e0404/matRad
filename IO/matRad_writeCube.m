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
%   datatype:                   MATLAB numeric datatype
%   meta:                       meta-information in struct. 
%                               Necessary fieldnames are:
%                               - resolution: [x y z]
%                               - datatype: numeric MATLAB-Datatype
%                               Optional:
%                               - axisPermutation (matRad default [2 1 3])
%                               - coordinateSystem (matRad default 'LPS')
%                               - imageOrigin (as used in DICOM)
%                               - dataUnit (i.e. Gy..)
%                               - dataName (i.e. dose, ED, ...)
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

%% Sanity checks
[filedir,filename,ext] = fileparts(filepath);
if ~exist(filedir,'dir')
    error(['Directory ' filedir ' does not exist!']);
end

%No Special characters in filename (allow only _ and alphanumeric
%characters
robustFileName = filename(isstrprop(filename,'alphanum') | filename == '_');
if ~strcmp(robustFileName,filename)
    warning(['Changing filename from ''' filename ''' to ''' robustFileName ''' to get rid of special characters!']);
    filepath = fullfile(filedir,[robustFileName ext]);
end

%% Prepare Metadata
%if the field is not set, we assume standard matRad x-y swap
if ~isfield(metadata,'axisPermutation')
    cube = permute(cube,[2 1 3]);
    metadata.axisPermutation = [1 2 3]; %Default Matlab axis permutation
end
%use the matrad coordinate system
if ~isfield(metadata,'coordinateSystem')
    metadata.coordinateSystem = 'LPS';  %Matlab coordinate system
end
%If there is no image origin set, center the image
imageExtent = metadata.resolution .* size(cube);
if ~isfield(metadata,'imageOrigin')
    metadata.imageOrigin = zeros(1,numel(imageExtent)) - (imageExtent/2);
end
%we can also store the center
if ~isfield(metadata,'imageCenter')
    metadata.imageCenter = metadata.imageOrigin + (imageExtent/2);
end


metadata.datatype = datatype;
    
%% Choose writer
%So far we only have an nrrd writer
switch ext
    case '.nrrd'
        matRad_writeNRRD(filepath,cube,metadata);
    case '.vtk'
        matRad_writeVTK(filepath,cube,metadata);
    case '.mha'
        matRad_writeMHA(filepath,cube,metadata);
    otherwise
        errordlg(['No writer found for extension "' ext '"']);
end

fprintf('File written to %s...\n',filepath);


end

