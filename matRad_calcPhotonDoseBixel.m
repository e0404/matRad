function dose = matRad_calcPhotonDoseBixel(SAD,m,betas,Interp_kernel1,...
                  Interp_kernel2,Interp_kernel3,radDepths,geoDists,...
                  latDistsX,latDistsZ)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad photon dose calculation for an individual bixel
% 
% call
%   dose = matRad_calcPhotonDoseBixel(SAD,Interp_kernel1,...
%                  Interp_kernel2,Interp_kernel3,radDepths,geoDists,...
%                  latDistsX,latDistsZ)
%
% input
%   SAD:                source to axis distance
%   m:                  absorption in water (part of the dose calc base
%                       data)
%   betas:              beta parameters for the parameterization of the 
%                       three depth dose components
%   Interp_kernel1/2/3: kernels for dose calculation
%   radDepths:          radiological depths
%   geoDists:           geometrical distance from virtual photon source
%   latDistsX:          lateral distance in X direction in BEV from central ray
%   latDistsZ:          lateral distance in Z direction in BEV from central ray
%
% output
%   dose:   photon dose at specified locations as linear vector
%
% References
%   [1] http://www.ncbi.nlm.nih.gov/pubmed/8497215
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define function_Di
func_Di = @(beta,x) beta/(beta-m) * (exp(-m*x) - exp(-beta*x)); 

% scale lateral distances to iso center plane
latDistsX = (latDistsX) ./ geoDists .* SAD;
latDistsZ = (latDistsZ) ./ geoDists .* SAD;
       
% Calulate lateral distances using grid interpolation.
lat1 = Interp_kernel1(latDistsX,latDistsZ);
lat2 = Interp_kernel2(latDistsX,latDistsZ);
lat3 = Interp_kernel3(latDistsX,latDistsZ);

% now add everything together (eq 19 w/o inv sq corr -> see below)
dose = lat1 .* func_Di(betas(1),radDepths) + ...
       lat2 .* func_Di(betas(2),radDepths) + ...
       lat3 .* func_Di(betas(3),radDepths);

% inverse square correction
dose = dose .* (SAD./geoDists(:)).^2;

% check if we have valid dose values
if sum(isnan(dose)) ||  sum(dose<0)>0
   error('Error in photon dose calc\n\n');
end
