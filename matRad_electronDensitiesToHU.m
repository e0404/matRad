function ct = matRad_electronDensitiesToHU(ct)
% matRad function to calculate recalculate HU values from equivalent 
% densities. This is done to provide downward compatability to previous
% matRad versions where HU values were not automatically saved during the
% import process. HU values can only be calculated if the HLUT is
% bijective.
%
% call
%   ct = matRad_electronDensitiesToHU(ct)
%
% input
%   ct: matRad ct struct containing cube and all additional information
%
% output
%   ct: ct struct with HU and equivalent density cube
%
% References
%   -
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



% load hlut
if isfield(ct,'hlut') % if hlut stored upon import use this one!
    hlut = ct.hlut;
else
    hlut = matRad_loadHLUT(ct);
end

% interpolate rel. electron dens. to HU based on lookup table
if isequal(hlut(:,2),unique(hlut(:,2))) && isequal(hlut(:,1),unique(hlut(:,1)))

    for i = 1:ct.numOfCtScen
        % Manual adjustments if ct data is corrupt. If some values are out of range
        % of the LUT, then these values are adjusted.
        if max(ct.cube{i}(:)) > max(hlut(:,2))
            warning('projecting out of range electron density values');
            ct.cube{i}(ct.cube{i}(:) > max(hlut(:,2))) = max(hlut(:,2));
        end
        if min(ct.cube{i}(:)) < min(hlut(:,2))
            warning('projecting out of range electron density values');
            ct.cube{i}(ct.cube{i}(:) < min(hlut(:,2))) = min(hlut(:,2));
        end

        ct.cubeHU{i} = interp1(hlut(:,2),hlut(:,1),ct.cube{i});
    end
    
else
    fprintf('Reconversion of HU values could not be done because HLUT is not bijective.\n');
end
