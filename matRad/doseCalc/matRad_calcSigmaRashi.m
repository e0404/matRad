function sigmaRashi = matRad_calcSigmaRashi(bdEntry, rangeShifter, SSD, radiationMode)
% calculation of additional beam broadening due to the use of range shifters
%
% call:
%   sigmaRashi = matRad_calcSigmaRashi(bdEntry,rangeShifter,SSD)
%   sigmaRashi = matRad_calcSigmaRashi(bdEntry,rangeShifter,SSD,radiationMode)
%
% input:
%   bdEntry:        base data entry for energy. Needs the field 'range'
%                   (ion range in water [mm]) or 'energy' (kinetic energy,
%                   [MeV] for protons / [MeV/u] for heavier ions)
%   rangeShifter:   structure defining range shifter geometry
%   SSD:            source to surface distance [mm]
%   radiationMode:  (optional) particle type, i.e., 'protons' (default),
%                   'helium', 'carbon', or 'oxygen'
%
% output:
%   sigmaRashi:   sigma of range shifter (to be added ^2) in mm
%
% The parametrization from [1,2] is fitted for protons. Heavier ions are
% handled by scaling to an equivalent proton with the same velocity
% profile: for an ion with charge z and mass number A, depths & ranges
% scale with z^2/A while the (Highland) scattering power scales with
% (z/A)^2, assuming pv scales with A at equal velocity.
%
% References
%   [1] https://www.ncbi.nlm.nih.gov/pubmed/12375823
%   [2] https://www.ncbi.nlm.nih.gov/pubmed/12701891
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    radiationMode = 'protons';
end

% particle charge and mass number for scaling to the proton-based
% parametrization
switch radiationMode
    case 'protons'
        z = 1;
        A = 1;
    case 'helium'
        z = 2;
        A = 4;
    case 'carbon'
        z = 6;
        A = 12;
    case 'oxygen'
        z = 8;
        A = 16;
    otherwise
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispError('Range shifter scattering not implemented for radiation mode ''%s''!', radiationMode);
end

% distance of range shifter to patient surface
rashiDist = SSD - rangeShifter.sourceRashiDistance;

% convert to mm
zS = rashiDist / 10; % [cm] (after division)
t  = rangeShifter.eqThickness / 10; % [cm] (after division)

% --> everything here is in [cm]

a1 = 0.21; % [1]
a2 = 0.77; % [1]
c0 = 0.0191027; % [cm]
c1 = 0.0204539; % [1]
alpha = 0.0022; % [cm MeV ^ (-p)] %
p = 1.77; % [1] % Exponent of range-energy relation

% (ion) range in water [cm]
if isfield(bdEntry, 'range')
    R = bdEntry.range / 10; % base data range is stored in [mm]
else % convert energy to range (Bragg-Kleeman with energy per nucleon)
    energy = bdEntry.energy;
    R = (A / z^2) * alpha * (energy^p);
end

% check if valid computation possible or range shifter to thick
if t / R >= 0.95
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Computation of range shifter sigma invalid. Range shifter is too thick.');
end

% relative depth in the range shifter (invariant under proton scaling)
s = t / R;

% range of the equivalent proton (same velocity profile as the ion)
Rp = (z^2 / A) * R;

% Improved HONG's Formula [1, Eq. 21]
sigmaT = (a1 * s + a2 * s^2) * (c0 + c1 * Rp);

% See [1, Eq. 11]
C1 = 13.6; % MeV
C2 = 0.038;

% See [1, Table 1 & Eq. 13]
LR_water   = 36.08;  % rad. length [cm] % cm
% rSPwater   = 1;     %rSP [rel. u.]
% PrSFwater = 1;      %lateral scaling factor, [rel.u.]
% rho_water = 1;      %density [g/cm^3]

LR_pmma = 40.55;     % rad. length [cm] % cm
rSP_pmma = 1.165;    % rSP [rel. u.]
rho_pmma = 1.19;     % density [g/cm^3]
PrSF_pmma = sqrt(rSP_pmma^3 * LR_water / LR_pmma * (1 + 2 * C2 * log(rSP_pmma * LR_water / LR_pmma))); % 1.1877;

LR_lexan = 41.46;    % rad. length [cm]
rSP_lexan = 1.165;    % rSP [rel. u.]
rho_lexan = 1.20;     % density [g/cm^3]

%{
L_air  = 30.4e3;    %rad. length [cm]
PSP_air = 1.07e-3;  %rSP [rel. u.]
rho_air = 1.205e-3; %density [g/cm^3]
PrSF_air = 990.81;  %lateral scaling factor, [rel.u.]
%}

F1part1 = (2 * C1^2 * alpha^(2 / p)) / (4 * LR_pmma * (2 / p - 1));
F1part2 = (((1 - s)^(2 - 2 / p) - 1) / (2 / p - 2) - s);
F1part3 = Rp^(2 - 2 / p);

F1 = F1part1 * F1part2 * F1part3;

F2part1 = (C1^2 * alpha^(2 / p)) / (4 * LR_pmma * (2 / p - 1));
F2part2 = ((1 - s)^(1 - 2 / p) - 1) / 1;
F2part3 = Rp^(1 - 2 / p);

F2 = F2part1 * F2part2 * F2part3;

% combine the moments, scaling each from the equivalent proton back to the
% ion (depth scaling z^2/A, scattering power scaling (z/A)^2)
sigmaProjSq = PrSF_pmma^2 * ((A / z^4) * sigmaT^2 + (1 / z^2) * F1 .* zS * rSP_pmma + (1 / A) * F2 * (zS * rSP_pmma).^2); % [cm ^ 2]

% <-- below in [mm]

sigmaRashi = 10 * sqrt(sigmaProjSq); % [mm] (after factor 10)
