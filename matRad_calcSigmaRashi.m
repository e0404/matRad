function sigmaRashi = matRad_calcSigmaRashi(energy,rangeShifter,SSD)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculation of additional beam broadening due to the use of range shifters (only for protons)
% 
% call
%   sigmaRashi = matRad_calcSigmaRashi(rangeShifter,SSD)
%
% input
%   energy:       initial particle energy
%   rangeShifter: structure defining range shifter geometry
%   SSD:          source to surface distance
%
% output
%   sigmaRashi:   sigma of range shifter (to be added ^2) in mm
%
% References
%   https://www.ncbi.nlm.nih.gov/pubmed/12375823
%   https://www.ncbi.nlm.nih.gov/pubmed/12701891
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

% convert energy to range
R = alpha * (energy ^ p);

% check if valid computation possible or range shifter to thick
if t / R >= 0.95
    error('Computation of range shifter sigma invalid.');
end

% Improved HONG's Formula
s = t / R;
sigmaT = (a1 * s + a2 * s ^ 2) * (c0 + c1 * R);

PSPpmma = 1.165;
PrSFpmma = 0.816;
% PSP_air = 0.00107;
% PrSF_air = 991;
C = 13.6; % MeV
L = 36.1; % cm

F1part1 = (2 * C ^ 2 * alpha ^ (2 / p)) / (4 * L * (2 / p - 1));
F1part2 = (((1 - s) ^ (2 - 2 / p) - 1) / (2 / p - 2) - s);
F1part3 = R ^ (2 - 2 / p);

F1 = F1part1 * F1part2 * F1part3;

F2part1 = (C ^ 2 * alpha ^ (2 / p)) / (4 * L * (2 / p - 1));
F2part2 = (((1 - s) ^ (1 - 2 / p) - 1)) / 1;
F2part3 = R ^ (1 - 2 / p);

F2 = F2part1 * F2part2 * F2part3;

sigmaProjSq = PrSFpmma ^ 2 * (sigmaT ^ 2 + F1 .* zS * PSPpmma + F2 * (zS * PSPpmma) .^ 2); % [cm ^ 2]

% <-- below in [mm]

sigmaRashi = 10 * sqrt(sigmaProjSq); % [mm] (after factor 10)

