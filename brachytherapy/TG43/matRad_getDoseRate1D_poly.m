function DoseRate = matRad_getDoseRate1D_poly(machine,r_mm)
% Calculation of radial dose Rate, interpolating using polynomes
%   1D dose rate formalism from Rivard et al. (2004): AAPM TG-43 update, 
%   page 639, Eq. 11:
%
% call 
%   DoseRate = matRad_getDoseRate1D_poly(machine,r_mm)
%
% input
%   machine:    TG43 information about the used seeds
%   r:          radial distance array, given in mm!
%
% output 
%   DoseRate:   size(r) array of dose Rate in cGy/h
%
% comment on dimensions / units
%   TG43 consensus data   cm, cGy, s
%   matRad                mm, Gy, s
%   output dimensions depend on the dimensions of air kerma strength
%   Sk, normallyi in cGy*cm^2/h)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2021 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% validate/ complete input arguments
if ~isfield(machine.data,'AnisotropyFactorRadialDistance')
    matRad_cfg.dispError('machine is missing field "AnisotropyFactorRadialDistance"...you might be trying to apply the TG43 2D formalism for basedata measured for 1D formalism') 
end
if ~isfield(machine.data,'AnisotropyFactorValue')
    matRad_cfg.dispError('machine is missing field "AnisotropyFactorValue"')
end
if ~isfield(machine.data,'lambda')
    matRad_cfg.dispError...
        ('machine is missing field "lambda" (dose constant in water)') 
end
if  machine.data.lambda < 0
    matRad_cfg.dispError('negative doseRate')
end
if min(r_mm,[],'all') < 0
    matRad_cfg.dispError('r contatins negative distances')
end
if ~isfield(machine.data,'SourceStrengthImplanted')
    machine.data.SourceStrengthImplanted = 1;
end 

%% arguments used during function
% r: radius (within this function all radii are given in cm)
r = 0.1*r_mm; 

% Sk: Air-kerma strength in  U...[1 U = 1 muGy*m^2/h) = 1 cGy*cm^2/h]
Sk = machine.data.SourceStrengthImplanted;

% lambda: Dose-rate constant in water (Lambda) in cGy/(h*U)
lambda = machine.data.lambda;

% L: length of line source in cm
L = machine.data.ActiveSourceLength;

% r0: reference radius in cm
r0 = 1;

% thet0: standard angle in degree
theta0 = 90;

% gLTab: Tabulated radial dose function \\ cell entry 1: radii; entry 2: values
% radii in cm, values in units of g(r0)
gLTab{1} = machine.data.RadialDoseDistance;
gLTab{2} = machine.data.RadialDoseValue;

% PhiAn: Tabulated anisotropy factor \\ cell entry 1: radii; entry 2: values
% radii in cm, values unitless
PhiAnTab{1} = machine.data.AnisotropyFactorRadialDistance;
PhiAnTab{2} = machine.data.AnisotropyFactorValue;

%% 1D formalism
% according to Rivard et al.: AAPM TG-43 update Eq. (11)
gL = matRad_radialDoseFunction(r,gLTab);
GL = matRad_geometryFunction(r,theta0,L);
GL0 = matRad_geometryFunction(r0,theta0,L);
PhiAn = matRad_anisotropyFactor1D(r,PhiAnTab, L);

DoseRate = Sk * lambda * GL./GL0 .* gL .* PhiAn;

end
