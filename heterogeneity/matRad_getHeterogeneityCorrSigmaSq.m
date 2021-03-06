function sigmaSq = matRad_getHeterogeneityCorrSigmaSq(WET,Pmod)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad calculation of Bragg peak degradation due to heterogeneities
% 
% call
%   sigmaSq = matRad_getHeterogeneityCorrSigmaSq(WET)
% 
% input
%   WET:        water equivalent thickness of heterogeneous structure, e.g. lung [mm]
%   Pmod:       modulation power [µm] (optional)
% 
% output
%   sigmaSq:    sigma squared [mm^2] for degradation of Bragg peak
% 
% References
%   -
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% a = 1.43 [mm/sqrt(cm)] from Riccardos thesis, fit MC well, 
% independent from proton energy, R80 stays at same position with degradation
% a = 1.60 [mm/sqrt(cm)] = 1.6/sqrt(10) [sqrt(mm)] fits the measurements better
% 
% Pmod = a^2; 150-750 micrometer for swine lung from Witt et al.
% Pmod = 256 [µm] is equivalent to a = 1.6 [mm/sqrt(cm)]
% Pmod = 204.5 [µm] to a = 1.43 [mm/sqrt(cm)]
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
%     Pmod = 256; % [µm]
    Pmod = 800;
%     fprintf('Default value for modulation power used: Pmod = 256 µm.\n');
end

sigmaSq = Pmod/1000 .* WET;
