function isoDoseContours = matRad_computeIsoDoseContours(doseCube,isoLevels)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function that computes all isodse contours (along each dose cube
% dimension)
%
% call
%   isoDoseContours = matRad_computeIsoDoseContours(doseCube,isolevels)
%
% input
%   doseCube            3D array containing the dose cube
%   isoLevels           iso dose levels (same units as doseCube)
%
% output
%   isoDoseContours     cell array containing the isolines along each dim
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

dim = size(doseCube);
isoDoseContours = cell(max(dim(:)),3);

minLevel = min(isoLevels(:));

for slice = 1:dim(1)
    if any(any(doseCube(slice,:,:) >= minLevel))
        isoDoseContours{slice,1} = contourc(squeeze(doseCube(slice,:,:)),isoLevels);
    end
end
for slice = 1:dim(2)
    if any(any(doseCube(:,slice,:) >= minLevel))
        isoDoseContours{slice,2} = contourc(squeeze(doseCube(:,slice,:)),isoLevels);
    end
end
for slice = 1:dim(3)
    if any(any(doseCube(:,:,slice) >= minLevel))
        isoDoseContours{slice,3} = contourc(squeeze(doseCube(:,:,slice)),isoLevels);
    end
end

end

