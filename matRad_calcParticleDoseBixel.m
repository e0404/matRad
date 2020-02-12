function dose = matRad_calcParticleDoseBixel(radDepths, radialDist_sq, sigmaIni_sq, baseData, stdWER)
% matRad visualization of two-dimensional dose distributions on ct including
% segmentation
% 
% call
%   dose = matRad_calcParticleDoseBixel(radDepths,radialDist_sq,SSD,focusIx,baseData)
%
% input
%   radDepths:      radiological depths
%   radialDist_sq:  squared radial distance in BEV from central ray
%   sigmaIni_sq:    initial Gaussian sigma^2 of beam at patient surface
%   baseData:       base data required for particle dose calculation
%
% output
%   dose:   particle dose at specified locations as linear vector
%
% References
%   [1] http://iopscience.iop.org/0031-9155/41/8/005
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

% add potential offset
depths = baseData.depths + baseData.offset;

% convert from MeV cm^2/g per primary to Gy mm^2 per 1e6 primaries
conversionFactor = 1.6021766208e-02;

if nargin < 5
   	stdCorr = false;
else
    stdCorr = true;
end


if ~isstruct(baseData.Z)
    if stdCorr
        for i = 1:size(stdWER,1)
            if stdWER(i) == 0
                kernel(i).values = 1;
            else
                tmp = normpdf(-3 * stdWER(i):0.05:3 * stdWER(i),0,stdWER(i));
                tmp = tmp / sum(tmp);
                kernel(i).values = tmp;
            end
        end
    end

    if ~isfield(baseData,'sigma')



        % interpolate depth dose, sigmas, and weights    
        if ~stdCorr
            X = matRad_interp1(depths,[conversionFactor*baseData.Z baseData.sigma1 baseData.weight baseData.sigma2],radDepths);
        else
            convBaseData = arrayfun(@(K) conv(conversionFactor*baseData.Z,K.values,'same'),kernel,'UniformOutput',false,'UniformOutput',false);
            X = [];
            for i = 1:size(radDepths)
                X = [X; matRad_interp1(depths,[convBaseData{i} baseData.sigma1 baseData.weight baseData.sigma2],radDepths(i))];
            end
        end


        % set dose for query > tabulated depth dose values to zero
        X(radDepths > max(depths),1) = 0;

        % compute lateral sigmas
        sigmaSq_Narr = X(:,2).^2 + sigmaIni_sq;
        sigmaSq_Bro  = X(:,4).^2 + sigmaIni_sq;

        % calculate lateral profile
        L_Narr =  exp( -radialDist_sq ./ (2*sigmaSq_Narr))./(2*pi*sigmaSq_Narr);
        L_Bro  =  exp( -radialDist_sq ./ (2*sigmaSq_Bro ))./(2*pi*sigmaSq_Bro );
        L = baseData.LatCutOff.CompFac * ((1-X(:,3)).*L_Narr + X(:,3).*L_Bro);

        dose = X(:,1).*L;
    else

        % interpolate depth dose and sigma
        X = matRad_interp1(depths,[conversionFactor*baseData.Z baseData.sigma],radDepths);

        %compute lateral sigma
        sigmaSq = X(:,2).^2 + sigmaIni_sq;

        % calculate dose
        dose = baseData.LatCutOff.CompFac * exp( -radialDist_sq ./ (2*sigmaSq)) .* X(:,1) ./(2*pi*sigmaSq);

     end
else
    % function handle for calculating lateral dose
    Gauss    = @(x,mu,SqSigma) 1./(sqrt(2*pi.*SqSigma)).*exp(-((x - mu).^2./(2.*SqSigma)));

    % function handle for calculating depth doses
    sumGauss = @(x,mu,SqSigma,w) (1./sqrt(2*pi*ones(numel(x),1) .* SqSigma') .* ...
        exp(-bsxfun(@minus,x,mu').^2 ./ (2* ones(numel(x),1) .* SqSigma' ))) * w;
    
    
    if ~isfield(baseData,'sigma')
        
        % interpolate depth dose, sigmas, and weights    
        %X = matRad_interp1(depths,[conversionFactor*baseData.Z.profileORG baseData.sigma1 baseData.weight baseData.sigma2],radDepths);
        X = matRad_interp1(depths,[baseData.sigma1 baseData.weight baseData.sigma2],radDepths);
        
        if ~stdCorr
            X = [conversionFactor * sumGauss(radDepths,baseData.Z.mean,baseData.Z.width.^2,baseData.Z.weight), X];
        else
            tmp = arrayfun(@(depth, std)  sumGauss(depth,baseData.Z.mean,baseData.Z.width.^2 + std.^2,baseData.Z.weight), radDepths, stdWER);
            X = [conversionFactor * tmp, X];
        end
        
        % set dose for query > tabulated depth dose values to zero
        X(radDepths > max(depths),1) = 0;

        % compute lateral sigmas
        sigmaSq_Narr = X(:,2).^2 + sigmaIni_sq;
        sigmaSq_Bro  = X(:,4).^2 + sigmaIni_sq;

        % calculate lateral profile
        L_Narr =  exp( -radialDist_sq ./ (2*sigmaSq_Narr))./(2*pi*sigmaSq_Narr);
        L_Bro  =  exp( -radialDist_sq ./ (2*sigmaSq_Bro ))./(2*pi*sigmaSq_Bro );
        L = baseData.LatCutOff.CompFac * ((1-X(:,3)).*L_Narr + X(:,3).*L_Bro);

        dose = X(:,1).*L;
    else

        % interpolate depth dose and sigma
        X = matRad_interp1(depths,[conversionFactor*baseData.Z baseData.sigma],radDepths);

        %compute lateral sigma
        sigmaSq = X(:,2).^2 + sigmaIni_sq;

        % calculate dose
        dose = baseData.LatCutOff.CompFac * exp( -radialDist_sq ./ (2*sigmaSq)) .* X(:,1) ./(2*pi*sigmaSq);

    end
end

% check if we have valid dose values
if any(isnan(dose)) || any(dose<0)
   error('Error in particle dose calculation.');
end 
