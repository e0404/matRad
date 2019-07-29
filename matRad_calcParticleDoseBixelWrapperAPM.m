function   [bixelDose] = matRad_calcParticleDoseBixelWrapperAPM(radDepths,latDistsX,latDistsZ,sigmaIni_sq,baseEntry,flagBioOpt,vTissueIndex,calcLET,currHeteroCorrDepths)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad particle pencil beam calculation wrapper
%
% call
%   dij = matRad_calcParticleDose(ct,stf,pln,cst)
%
% input
%   radDepths:        radiological depths
%   latDistsX:        lateral distance to the central beam axis in x-direction
%   latDistsZ:        lateral distance to the central beam axis in z-direction
%   sigmaIni_sq:      initial Gaussian sigma^2 of beam at patient surface
%   baseEntry:        base data required for particle dose calculation
%   flagBioOpt        boolean determining biological or physical optimization
%
% output
%   bixelDose:        struct holding the deposited dose in x,y,z as well as
%                     the total dose
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% function handle for calculating depth doses
sumGauss = @(x,mu,SqSigma,w) ((1./sqrt(2*pi*ones(numel(x),1) * SqSigma') .* ...
    exp(-bsxfun(@minus,x,mu').^2 ./ (2* ones(numel(x),1) * SqSigma' ))) * w);

if ~exist('calcLET','var')
    calcLET = false;
end

if exist('currHeteroCorrDepths','var')
    calcAPM = true;
end

% perform bixel dose calculation depending on the given base data format
if calcAPM
    [bixelDose] = matRad_calcParticleDoseBixelAPM(radDepths,latDistsX,latDistsZ,sigmaIni_sq,baseEntry,flagBioOpt,vTissueIndex,currHeteroCorrDepths);
else
    rad_distancesSq    = latDistsX.^2 + latDistsZ.^2;
    bixelDose.physDose = matRad_calcParticleDoseBixel(radDepths,rad_distancesSq,sigmaIni_sq,baseEntry);
end

if calcLET
    depths = baseEntry.depths + baseEntry.offset;
    if ~isstruct(baseEntry.LET)
        bixelDose.LET_ij = matRad_interp1(depths,baseEntry.LET,radDepths);
    else
        bixelDose.LET_ij   = sumGauss(radDepths, baseEntry.LET.mean, (baseEntry.LET.width).^2,baseEntry.LET.weight);
    end
    bixelDose.LET_ij(isnan(bixelDose.LET_ij)) = 0;
end

end