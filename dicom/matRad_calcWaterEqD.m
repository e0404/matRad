function ct = matRad_calcWaterEqD(ct, radiationMode)
% matRad function to calculate the equivalent densities from a dicom ct 
% that originally uses intensity values
%
% call
%   ct = matRad_calcWaterEqD(ct, radiationMode)
%
% input
%   ct:             ct containing a cubeHU to compute rED/rSP values from
%   radiationMode:  radiationMode as character array (e.g. 'photons') since matRad 3.
%                   Can also be a pln-struct for downwards compatibility
%
% output
%   ct: ct struct with cube with relative _electron_ densities stored in
%   ct.cube
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


matRad_cfg = MatRad_Config.instance();

% load hlut
hlut = matRad_loadHLUT(ct, radiationMode);  

for i = 1:ct.numOfCtScen

    % Manual adjustments if ct data is corrupt. If some values are out of range
    % of the LUT, then these values are adjusted.
    if max(ct.cubeHU{i}(:)) > max(hlut(:,1))
        matRad_cfg.dispWarning('projecting out of range HU values');
        ct.cubeHU{i}(ct.cubeHU{i} > max(hlut(:,1))) = max(hlut(:,1));
    end
    if min(ct.cubeHU{i}(:)) < min(hlut(:,1))
        matRad_cfg.dispWarning('projecting out of range HU values');
        ct.cubeHU{i}(ct.cubeHU{i} < min(hlut(:,1))) = min(hlut(:,1));
    end

    % interpolate HU to relative electron density or relative stopping power based on lookup table
    ct.cube{i} = reshape(matRad_interp1(hlut(:,1),hlut(:,2),double(ct.cubeHU{i}(:))),ct.cubeDim);

end

% save hlut
ct.hlut = hlut;

end
