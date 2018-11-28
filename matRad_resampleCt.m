
function ctNew = matRad_resampleCt(ct, res)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to resample cube down or up 
% 
% call
%   ctNew = matRad_resampleCt(ct, res)
%
% input
%   ct:             matRad ct structure
%   res:            new resolution of the cubes(mm)
%
% output
%   ct:             matRad ct struct. Note that this 3D matlab array 
%                   contains water euqivalen electron denisities.
%                   Hounsfield units are converted using a standard lookup
%                   table in matRad_calcWaterEqD
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
ctNew.resolution = res ;
ctNew.x = ct.x(1):ctNew.resolution(1):ct.x(end);
ctNew.y = ct.y(1):ctNew.resolution(2):ct.y(end);
ctNew.z = ct.z(1):ctNew.resolution(3):ct.z(end);
ctNew.cubeDim = [numel(ctNew.x) numel(ctNew.y) numel(ctNew.z)];
ctNew.numOfCtScen = ct.numOfCtScen;
ctNew.dicomInfo = ct.dicomInfo;
ctNew.dicomMeta = ct.dicomMeta;
ctNew.timeStamp = datestr(now);
ctNew.hlut = ct.hlut;

[X,Y,Z] = meshgrid(ct.x,ct.y,ct.z);
[Xq,Yq, Zq] = meshgrid(ctNew.x, ctNew.y, ctNew.z);
% ct cube
for i = 1:numel(ct.cube)
ctNew.cube{i} = interp3(X,Y,Z, ct.cube{i}, Xq,Yq,Zq);
% ct HU
ctNew.cubeHU{i} = interp3(X,Y,Z, ct.cubeHU{i}, Xq,Yq,Zq);

end
end
