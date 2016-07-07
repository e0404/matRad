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
%   meta:                       meta-information in struct. Fieldnames are
%                               - datatype
%                               - resolution [x y z]
%                               - coordinates [x y z] cube center
%                               - orientation (default ..)
%
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

%% x-y permutation from matRad coordinate system
%cube = permute(cube,[2 1 3]);

%% Prepare Metadata
if ~isfield(metadata,'axisPermutation')
    metadata.axisPermutation = [2 1 3]; %Default Matlab axis permutation
end
if ~isfield(metadata,'coordinateSystem')
    metadata.coordinateSystem = 'LPS';  %Matlab coordinate system
end

%% Choose writer
switch ext
    case '.nrrd'
        matRad_writeNRRD(filepath,cube,datatype,metadata);
    otherwise
        errordlg(['No writer found for extension "' ext '"']);
end


end

