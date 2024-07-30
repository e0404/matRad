function intDose = matRad_calcIntEnergy(dose,ct,pln)
% matRad function to compute the integral energy in MeV for a dose cube
%
% call
%   intDose = matRad_calcIntEnergy(dose,ct,pln)
%
% input
%   dose:   3D matlab array with dose e.g. resultGUI.physicalDose
%   ct:     matRad ct struct
%   pln:    matRad pln struct
%
% output
%   intDose: integral dose in MeV
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
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

% make HU to RSP/electron density conversion if not available
if ~isfield(ct,'cube')
    ct = matRad_calcWaterEqD(ct, pln);
end
    
% conversion factor for J to eV
electronCharge = 1.60217662e-19; % [C]
convFac = 1/electronCharge;

% hints: ct.cube corresponds to g/cm^3 --> /1000 to have kg/cm^3
%        ct.resolution given in mm --> /1000 to have cm
%        conversion to MeV from eV --> /1000000
intDose = convFac * sum(dose(:).*ct.cube{1}(:)) ...
         * ct.resolution.x*ct.resolution.y*ct.resolution.z/1e12;

if nargout == 0
    matRad_cfg.dispInfo('Integral energy in dose cube = %.4g MeV\n',intDose);
end

