function dose = matRad_calcParticleDoseBixel(radDepths,radialDist_sq,SSD,focusIx,baseData)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad visualization of two-dimensional dose distributions on ct including
% segmentation
% 
% call
%   dose = matRad_calcParticleDoseBixel(radDepths,radialDist_sq,SSD,focusIx,baseData)
%
% input
%   radDepths:      radiological depths
%   radialDist_sq:  squared radial distance in BEV from central ray
%   SSD:            source to surface distance
%   focusIx:        index of focus to be used
%   baseData:       base data required for particle dose calculation
%
% output
%   dose:   particle dose at specified locations as linear vector
%
% References
%   [1] http://iopscience.iop.org/0031-9155/41/8/005
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

% range shift
depths = baseData.depths + baseData.offset;

% convert from MeV cm^2/g per primary to Gy mm^2 per 1e6 primaries
conversionFactor = 1.6021766208e-02;

 % calculate initial focus sigma
SigmaIni = interp1(baseData.initFocus(focusIx).dist,baseData.initFocus(focusIx).sigma,SSD);

if ~isfield(baseData,'sigma')
    
    % interpolate depth dose, sigmas, and weights    
    X = interp1(depths,[conversionFactor*baseData.Z baseData.sigma1 baseData.weight baseData.sigma2],radDepths,'linear');
    
    % compute lateral sigmas
    sigmaSq_Narr = X(:,2).^2 + SigmaIni^2;
    sigmaSq_Bro  = X(:,4).^2 + SigmaIni^2;
    
    % calculate lateral profile
    L_Narr =  exp( -radialDist_sq ./ (2*sigmaSq_Narr))./(2*pi*sigmaSq_Narr);
    L_Bro  =  exp( -radialDist_sq ./ (2*sigmaSq_Bro ))./(2*pi*sigmaSq_Bro );
    L = baseData.LatCutOff.CompFac * ((1-(X(:,3))).*L_Narr) + (X(:,3).*L_Bro);

    dose = X(:,1).*L;
else
    
    % interpolate depth dose and sigma
    X = interp1(depths,[conversionFactor*baseData.Z baseData.sigma],radDepths,'linear');

    %compute lateral sigma
    sigmaSq = X(:,2).^2 + SigmaIni^2;
    
    % calculate dose
    dose = baseData.LatCutOff.CompFac * exp( -radialDist_sq ./ (2*sigmaSq)) .* X(:,1) ./(2*pi*sigmaSq);
    
 end
 
% check if we have valid dose values
if any(isnan(dose)) || any(dose<0)
   error('Error in particle dose calculation.');
end 
