function [cube, metadata] = matRad_readCube(filename)
% matRad Cube read wrapper
% determines the extension and assigns the appropriate reader to it
% 
% call
%   matRad_readCube(filename)
%
% input
%   filename:   full path of the file
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

if ~exist(filename,'file')
    error(['File ' filename ' does not exist!']);
end

[pathstr,name,ext] = fileparts(filename);

switch ext
    case {'.nrrd','.NRRD'}
        disp(['Reading NRRD: ' filename '...']);
        [cube, metadata] = matRad_readNRRD(filename);
        disp('Done!');
    otherwise
        error(['Extension ' ext ' not (yet) supported!']);
end
metadata.name = name;
metadata.path = pathstr;

end

